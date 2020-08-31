
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

Y = wrandnp(Float32, 2,5)
@test eltype(Y.array) == eltype(weights(Y)) == Float32

Z = wrandn(big, 2,5)
@test eltype(Z.array) == eltype(weights(Z)) == BigFloat

pca = wPCA(R)
@test size(pca(sobol(3,10))) == (2,10) # pca can be applied

sPCA(R)
@test size(rPCA(sobol(3,20))) == (2,20) # was saved

using SliceMap

@test mapslices(x -> 2 .* x, S).array ≈ 2 .* S.array

@test mapcols(x -> 2 .* x, S).array ≈ 2 .* S.array
@test MapCols(x -> 2 .* x, S).array ≈ 2 .* S.array

mktempdir() do path
    w1 = Weighted([1 2; 3 4], [0.5, 1.0], WeightOpt(aname="a", wname="w"))

    f1csv = joinpath(path, "w1.csv")
    save(w1, f1csv)
    w1csv = load(f1csv)

    @test w1csv.array == w1.array
    @test w1csv.weights == weights(w1) # normalised
    @test w1csv.opt.aname == "w1" # not saved

    f1json = joinpath(path, "w1.json")
    save(w1, f1json)
    w1json = load(f1json)

    @test w1json.array == w1.array
    @test w1json.weights == w1.weights
    @test w1json.opt.aname == "a"
end
