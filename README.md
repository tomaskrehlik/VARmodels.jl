# VARmodels

[![Build Status](https://travis-ci.org/tomaskrehlik/VARmodels.jl.svg?branch=master)](https://travis-ci.org/tomaskrehlik/VARmodels.jl)

This package allows to fit various autoregressive models. So far it only knows simple vector autoregressions (VAR) and restricted VAR using EGLS method. The package is under development and by the end of September it should include also vector error correction models (VECM). The VECM models will depend on a separate unit-root testing package.

## Starting example

First, read in some data, I am using for example Realized Volatilities from [here](http://realized.oxford-man.ox.ac.uk/data/download) with removed headers.

````julia
data = readcsv("data.csv")

# Only use three of them
e = varEstimate(data[:,[ 2,30, 80]], 4, "Const")

# The coefficients
e.C

# The impuls responses
Psi(e, 10)
using PyPlot
plot([e.Psi[1,1,i] for i=1:10])

# Restriction of coefficients to zero
# Matrix that has the same size as the coefficient matrix
# Put false for which coeffcient you do not want to appear
rest = fill(true, 3, 13)
rest[1, 5] = false

# Estimate, restriction matrix, if you want to fit using EGLS
restrictVAR2(e, rest, true)
````

## TO-DO

- Finish generalised restrictions. (See Lutkepohl, H. 2007, pg. 194) Mostly come up with function that generates suitable matrix `R` from output of `@restrictions` macro. So far undocumented macro, but see examples in the source code.
- Add standard errors, even though it is not really neccessary for VAR, as usually the impulse response std errors are important
- Add bootstrapped errors for the impulse responses. 
- Add tests for anything. Basically just rewrite formulas using the `λ` function.
- Start doing VECM, focusing first on unit-root tests.
