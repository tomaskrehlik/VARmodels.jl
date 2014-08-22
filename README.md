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

````