##################################################################
#                                             twalkexample.jl
#
#  Example for the twalk implementation in Julia,
#  Created by J Andres Christen, jac at cimat.mx .
#  Julia implementation by Nicolás Kuschinski
#
#  See http://www.cimat.mx/~jac/twalk/ for more details.
#  and for the Python, R, Matlab and C implementations.
#
##################################################################

### Differences between jtwalk and pytwalk versions:
## The Julia version does not perform output analysis
## Does not estimate the remaining time
## and does not include specialized plotting functions.

## Plotting is easy in julia! The following example
## includes sample plots using the Plots.jl plotting library


using JTwalk


####################################################################
######### Product of Independent Exponentilas #########

###   In the jtwalk object, we define the objective function with an
###   energy function U.
###   The support is defined in a separate function, and is defined in
###   jtwalk as Supp.

###   NOTE: The objective function is defined as -log of the target
###   density function

lambdas = [ 1., 2., 3., 4., 5.]

function ExpU(x)
    # -log of a product of exponentials
    return sum(x .* lambdas)
end

function ExpSupp(x)
	return all(0. .< x)
end


###   The following initializes the TWalk object
###   The dimension of the parameter space is n
Exp = jtwalk( n=5, U=ExpU, Supp=ExpSupp)

#### This runs the T-Walk

Run!(Exp, T=50000, x0=30*ones(5), xp0=40*ones(5))

##### The trajectory of the chain is stored into Exp.Output & Exp.Output_u
##### These are (n+1) by T matrices, with one row per iteration.
##### Each row has the output followed by the energy function evaluated
##### at the corresponding location. Both matrices are valid MCMC chains
##### but THEY ARE NOT INDEPENDENT.

### The following shows how to use the Plots.jl plotting library
### to make relevant plots. In julia even plotting complex thinning can
### be done with simple one line instructions.

using Plots

### The following produces a trace plot:

plot(Exp.Output[:,end])

### The following shows a scatter plot of the first and third variables with
### thinning of 1000 iteration burn-in and 250 iteration autocorrelation time:

scatter(Exp.Output[1000:250:end,1], Exp.Output[1000:250:end,3])

### The following shows a histogram of the 4th variable with the same thinning:

histogram(Exp.Output[1000:250:end,4])
