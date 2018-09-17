

This package defines a WeightedMatrix, a struct with a matrix and a vector of weights corresponding to its columns.
I'm in the process of separated this out & upgrading to Julia 0.7.

What's working:
* sobol, wrand, etc.
* unique (with right GroupSlices version)
* hcat
* PCA (but with warnings from MultivariateStats)
* plots
* saving

Still to do:
* slicing
* symm
* broadcasting

[![Build Status](https://travis-ci.org/mcabbott/WeightedArrays.jl.svg?branch=master)](https://travis-ci.org/mcabbott/WeightedArrays.jl)
