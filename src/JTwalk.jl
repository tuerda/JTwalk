################################################
#           jtwalk
#   Twalk in Julia.
#   Twalk by Andres Christen & Colin Fox
#   Julia implementation by Nicolas Kuschinski
#
################################################

#   This is able to run the twalk algorithm, but:
#     -does not estimate remaining time
#     -does not perform output analysis
#       - does not include a RWMH.
#   Run and onemove methods are exported.
#   And can be accessed independently.

# This is mostly just translating the python version of twalk,
# some optimization has been done for Julia:
# Type instability has been resolved everywhere except in Energy! SetupInitialValues! and onemove.
# SetupInitialValues! and Energy! do not matter because they are only used once.

# Type instability in onemove has been isolated to the invocation of the Supp and U functions.
# Profiling and analysis done by José Zubieta indicates that this type instability is not serious
# and that the performance hit is minimal





###############################################
# NOTE:
# This uses nothing outside of basic julia
# therefore no external libraries are required.
###############################################



module JTwalk
    export Twalk
    export jtwalk
    export Run!
    export onemove


#Auxilliary function:

    function SqrNorm(x::Array{Float64,1})
        sum(x .* x)
    end

    #Twalk structure
    mutable struct Twalk
        n::Int64
        U::Function
        Supp::Function
        t::Int64
        u::Float64
        up::Float64
        aw::Float64
        at::Float64
        pphi::Float64
        Output::Array{Float64,2}
        Output_u::Array{Float64,2}
        T::Int64
        Fw::Array{Float64,1}
        LikelihoodEnergy::Function
        PriorEnergy::Function
        ll_e::Float64
        prior_e::Float64
        nphi::Int64
        phi::Array{Bool,1}
        sigma::Float64
        x::Array{Float64,1}
        xp::Array{Float64,1}
    end

    function Energy!(T::Twalk, x::Array{Float64,1})
        T.ll_e=T.LikelihoodEnergy(x)
        T.prior_e=T.PriorEnergy(x)
        T.t*T.ll_e+T.prior_e
    end

    #Initializer All arguments are keywords
    function jtwalk(;n::Int64, U=x->x[1]^2/2, Supp=x->true, t=-1, u=0, up=0, w=x->0., ww=[0., 0.4918, 0.4918, 0.0082, 0.0082], aw=1.5, at=6., n1phi=4.)
        if t>0  #Penalized likelihood
            LikelihoodEnergy=u
            PriorEnergy=w
            selft=t
        else
            PriorEnergy=x->0
            LikelihoodEnergy=U
            selft=1
        end
        pphi=min(n,n1phi)/n ##Prob of choosing each parameter
        # 
        #Create the object first & then modify it a little. Matches python version more closely.
        twalkinitial=Twalk(n, U, Supp, t, 0., up, aw, at, pphi, zeros(2,2), zeros(2,2), 1, cumsum(ww), LikelihoodEnergy, PriorEnergy, 0., 0., 0,[true],0.,[0],[0])
        # 
        twalkinitial.U=x->Energy!(twalkinitial,x)
        # 
        twalkinitial
    end

    function SetupInitialValues!(twalk::Twalk, x0::Array{Float64,1}, xp0::Array{Float64,1})
        #Check if x0 & xp0 are in the support
        (any(abs.(x0-xp0) .<=0)) && error("Not all entires of the initial values are different")
        (twalk.Supp(x0)) || error("x0 out of support")
        (twalk.Supp(xp0)) || error("xp0 out of support")
        u=twalk.U(x0)
        up=twalk.U(xp0)
        [true u up]
    end

    ###############################
    #Run the twalk

    function Run!(twalk::Twalk;T::Int64,x0::Array{Float64,1},xp0::Array{Float64,1},t=1)
        #Run the twalk. T=number of iterations, x0 & xp0 are start points.
        twalk.t=t
        # 
        (rt, u, up)=SetupInitialValues!(twalk,x0,xp0)
        Bool(rt) || return 0
        x=x0 #last values used
        xp=xp0
        twalk.Output = zeros(T+1,twalk.n+1)
        twalk.Output_u=zeros(T+1,twalk.n+1)
        twalk.T = T+1
        kercall = zeros(6)
        # 
        twalk.Output[1, 1:twalk.n] = copy(x)
        twalk.Output[1, twalk.n+1]=u
        twalk.Output_u[1]=twalk.ll_e
        #
        j1=1
        j=0
        # 
        #Sampling
        for it in 1:T
            y, yp, ke, A, u_prop, up_prop = onemove(twalk, x, u, xp, up)
            kercall[ke]+=1
            kercall[6]+=1
            if rand()<A
                x=copy(y)
                u=u_prop
                xp=copy(yp)
                up=up_prop
            end 
            #
            #To retrieve current values
            twalk.x=x
            twalk.xp=xp
            twalk.u=u
            twalk.up=up
            # 
            twalk.Output[it+1,1:twalk.n]=copy(x)
            twalk.Output[it+1,twalk.n+1]=u
            twalk.Output_u[it+1]=twalk.ll_e
        end
        # 
        return 1
    end

    ###############################
    #One move of twalk

    function onemove(twalk::Twalk, x::Array{Float64,1}, u::Float64, xp::Array{Float64,1}, up::Float64)
        #One move of the twalk
        n=twalk.n
        U=twalk.U
        Supp=twalk.Supp
        Fw=twalk.Fw
        ker=rand()
        ke=1
        u_prop::Float64=0.0
        up_prop::Float64=0.0
        A::Float64=0.0
        # 
        #kernel nothing exchange x with xp, ~ used
        if ((0<=ker) & (ker<Fw[1]))
            ke=0
            y=copy(xp)
            up_prop=u
            yp=copy(x)
            u_prop=up
            A=1.0
        end
        # 
        #The walk move
        if (Fw[1] <= ker < Fw[2])
            ke=1
            dir=rand()
            if(dir<0.5) #x as pivot
                yp=SimWalk(twalk,xp,x)
                y=copy(x)
                u_prop=u
                if (Supp(yp)& all(yp.!=y))
                    up_prop=U(yp)
                    A=exp(up-up_prop)
                else
                    up_prop=0.0
                    A=0.0
                end
            else #xp as pivot
                y=SimWalk(twalk,x,xp)
                yp=copy(xp)
                up_prop=up
                if(Supp(y) & (all(yp.!=y)))
                    u_prop=U(y)
                    A=exp(u-u_prop)
                else
                    u_prop=0.0
                    A=0.0
                end
            end
        end
        # 
        #The traverse move
        if (Fw[2] <= ker < Fw[3])
            ke=2
            dir=rand()
            if(dir<0.5) #x as pivot
                beta=Simbeta(twalk)
                yp=SimTraverse(twalk,xp,x,beta)
                y=copy(x)
                u_prop=u
                if Supp(yp)
                    up_prop=U(yp)
                    (twalk.nphi==0) ? A=1.0 : A=(exp((up-up_prop) + (twalk.nphi-2)*log(beta)))
                else
                    up_prop=0.0
                    A=0.0
                end
            else #xp as pivot
                beta=Simbeta(twalk)
                y=SimTraverse(twalk,x,xp,beta)
                yp=copy(xp)
                up_prop=up
                if Supp(y)
                    u_prop=U(y)
                    (twalk.nphi==0) ? A=1.0 : A=(exp((u-u_prop) + (twalk.nphi-2)*log(beta)))
                else
                    u_prop=0.0
                    A=0.0
                end
            end
        end
        # 
        #The blow move
        if (Fw[3] <= ker < Fw[4])
            ke=3
            dir=rand()
            if(dir<0.5) #x as pivot
                yp=SimBlow(twalk,xp,x)
                y=copy(x)
                u_prop=u
                if (Supp(yp)& all(yp.!=x))
                    up_prop=U(yp)
                    W1=GBlowU(twalk,yp,xp,x)
                    W2=GBlowU(twalk,xp,yp,x)
                    A=exp((up-up_prop)+(W1-W2))
                else
                    up_prop=0.0
                    A=0.0
                end
            else #xp as pivot
                y=SimBlow(twalk,x,xp)
                yp=copy(xp)
                up_prop=up
                if(Supp(y) & (all(y.!=xp)))
                    u_prop=U(y)
                    W1=GBlowU(twalk,y,x,xp)
                    W2=GBlowU(twalk,x,y,xp)
                    A=exp((u-u_prop)+(W1-W2))
                else
                    u_prop=0.0
                    A=0.0
                end
            end
        end
        # 
        #The Hop move
        if (Fw[4] <= ker < Fw[5])
            ke=4
            dir=rand()
            if(dir<0.5) #x as pivot
                yp=SimHop(twalk,xp,x)
                y=copy(x)
                u_prop=u
                if (Supp(yp)& all(yp.!=x))
                    up_prop=U(yp)
                    W1=GHopU(twalk,yp,xp,x)
                    W2=GHopU(twalk,xp,yp,x)
                    A=exp((up-up_prop)+(W1-W2))
                else
                    up_prop=0.0
                    A=0.0
                end
            else #xp as pivot
                y=SimHop(twalk,x,xp)
                yp=copy(xp)
                up_prop=up
                if(Supp(y) & (all(y.!=xp)))
                    u_prop=U(y)
                    W1=GHopU(twalk,y,x,xp)
                    W2=GHopU(twalk,x,y,xp)
                    A=exp((u-u_prop)+(W1-W2))
                else
                    u_prop=0.0
                    A=0.0
                end
            end
        end
        return [y,yp,ke,A,u_prop,up_prop]
    end

    ###############################
    #Auxiliary functions for the kernels

    #used by the walk kernel
    function SimWalk(twalk::Twalk,x::Array{Float64,1},xp::Array{Float64,1})
        aw=twalk.aw
        n=twalk.n
        phi=(rand(n).<twalk.pphi)
        twalk.nphi=sum(phi)
        z=zeros(n)
        for i in 1:n
            if phi[i]
                u=rand()
                z[i]=(aw/(1+aw))*(aw*u^2+2*u-1)
            end
        end
        return x+((x-xp) .* z)
    end

    #used by the traverse kernel
    function Simbeta(twalk::Twalk)
        at=twalk.at
        if rand()<((at+1)/(2*at))
            return exp(1/(at+1)*log(rand()))
        else
            return exp(1/(1-at)*log(rand()))
        end
    end

    function SimTraverse(twalk, x::Array{Float64,1}, xp::Array{Float64,1}, beta::Float64)
        n=twalk.n
        phi=(rand(n).<twalk.pphi)
        twalk.nphi=sum(phi)
        rt=copy(x)
        for i in 1:n
            if phi[i]
                rt[i]=xp[i]+beta*(xp[i]-x[i])
            end
        end
        return rt
    end

    #used by the blow kernel
    function SimBlow(twalk::Twalk,x::Array{Float64,1},xp::Array{Float64,1})
        n=twalk.n
        twalk.phi=rand(n).<twalk.pphi
        twalk.nphi=sum(twalk.phi)
        twalk.sigma=maximum(twalk.phi .* abs.(xp-x))
        rt=copy(x)
        for i in 1:n
            if twalk.phi[i]
                rt[i]=xp[i]+twalk.sigma*randn()
            end
        end
        return rt
    end

    function GBlowU(twalk::Twalk,h,x::Array{Float64,1},xp::Array{Float64,1})
        nphi=twalk.nphi
        twalk.sigma=maximum(twalk.phi .* abs.(xp-x))
        if (nphi>0)
            return(nphi/2.0)*log(2*pi) + nphi*log(twalk.sigma)+0.5*SqrNorm(h-xp)/(twalk.sigma^2)
        else
            return 0.0
        end
    end

    #used by the hop kernel
    function SimHop(twalk::Twalk,x::Array{Float64,1},xp::Array{Float64,1})
        n=twalk.n
        twalk.phi=rand(n).<twalk.pphi
        twalk.nphi=sum(twalk.phi)
        twalk.sigma=maximum(twalk.phi .* abs.(xp-x))/3
        rt=copy(x)
        for i in 1:n
            if twalk.phi[i]
                rt[i]=x[i]+twalk.sigma*randn()
            end
        end
        return rt
    end

    function GHopU(twalk::Twalk,h,x::Array{Float64,1},xp::Array{Float64,1})
        nphi=twalk.nphi
        twalk.sigma=maximum(twalk.phi .* abs.(xp-x))/3
        if (nphi>0)
            return(nphi/2.0)*log(2*pi) + nphi*log(twalk.sigma)+0.5*SqrNorm(h-x)/(twalk.sigma^2)
        else
            return 0.0
        end
    end

end #module
