
using WeightedArrays
using Test

R = Weighted(rand(3,10), rand(10))

@test size(R) == (3,10)
@test lastlength(R) == 10

@test sum(R.weights) ≈ 1 # normalised in construction

R.weights .+= 1
@test sum(weights(R)) ≈ 1 # normalised on reading

normalise!(R)
@test sum(R.weights) ≈ 1 # restored

@test tanh(R).opt.aname == "tanh-θ"
@test tanh(R).opt.lo == -1

S = sobol(2,9)

@test log(S).opt.clamp == false

@test size(unique(S + S)) == size(S) # hcat & unique

@test (2S).weights[1] ≈ 2/9
@test (2 .* S).array ≈ 2 .* (S.array)

X = xgrid(2,0:0.1:1)
X = near(X, S, 0.2)
@test size(X) == (2,78) # near

pca = wPCA(R)
@test size(pca(sobol(3,10))) == (2,10) # pca can be applied

sPCA(R)
@test size(rPCA(sobol(3,20))) == (2,20) # was saved

using SliceMap

@test mapslices(x -> 2 .* x, S).array ≈ 2 .* S.array

@test mapcols(x -> 2 .* x, S).array ≈ 2 .* S.array
@test MapCols(x -> 2 .* x, S).array ≈ 2 .* S.array
