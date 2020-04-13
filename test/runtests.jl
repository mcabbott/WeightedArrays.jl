
using WeightedArrays
using Test

R = Weighted(rand(3,10), rand(10))

@test size(R) == (3,10)
@test lastlength(R) == 10

@test sum(R.weights) ≈ 1 # normalised in construction

R.weights .+= 1
sum(weights(R)) == 1 # normalised on reading

normalise!(R)
@test sum(R.weights) ≈ 1 # restored

S = sobol(2,9)

@test size(unique(S + S)) == size(S) # hcat & unique

X = xgrid(2,0:0.1:1)
X = near(X, S, 0.2)
@test size(X) == (2,78) # near


pca = sPCA(R)
@test size(pca(sobol(3,10))) == (2,10) # pca can be applied

@test size(rPCA(sobol(3,20))) == (2,20) # was saved
