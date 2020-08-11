
export Weighted, WeightedMatrix, ClampedWeighted, UnClampedWeighted, WeightOpt
export array, weights, trace, lastlength, aname, wname, unclamp,
    normalise, normalise!, unnormalise, unnormalise!,
	totweight, wcopy, wcopy!, wglue, flatten, flatcopy, flatcopy!, svecs, CatView
export AbsVec, AbsMat, AbsArray

Base.@kwdef struct WeightOpt
    norm::Bool = true
    clamp::Bool = false
    lo::Float64 = 0.0
    hi::Float64 = 1.0
    aname::String = "θ"
    wname::String = "p(θ)"
    like::Bool = false
    trace = Float64[]
end

struct Weighted{T1,T2}
    array::T1
    weights::T2
    opt::WeightOpt

    function Weighted(array::T1, weights::T2, opt::WeightOpt) where {T1<:AbstractArray,T2<:AbstractVector}
        size(array,ndims(array)) == length(weights) || error("length of weights much match size of last dimension of array")
        new{T1, T2}(array, weights, opt)
    end
end

const WeightedArray = Weighted
const WeightedVector = Weighted{<:AbstractVector}
const WeightedMatrix = Weighted{<:AbstractMatrix}

const MaybeWeightedArray = Union{AbstractArray, WeightedArray}
const MaybeWeightedVector = Union{AbstractVector, WeightedVector}
const MaybeWeightedMatrix = Union{AbstractMatrix, WeightedMatrix}

const AbsVec = AbstractVector
const AbsMat = AbstractMatrix
const AbsArray = AbstractArray

const def_opt = WeightOpt()

##### Constructors

using LinearAlgebra

"""
    x = Weighted(array, weights)

Created a Weighted Array, with given weights associated to the last dimension of array,
e.g. to columns of a matrix.

Keyword `norm=true` specifies that `weights(x)` sum to 1 always;
this is enforced on construction. Calling `normalise!(x)` will re-enforce this,
e.g. after scaling `2x` or `hcat()`-ing `x + y` or mutating.

    x = Weighted(array, weights, lo, hi)

Here every `θ ∈ array` is constrained `lo ≦ θ ≦ hi`, the same for all dimensions.
(Weights are always `0 ≦ λ < ∞`.)

Calling `clamp!(x)` will re-enforce these box constraints,
for instance after mutating `x[1,2] *= 3`.

    x = Weighted(array, weights, opt)

This constructor with `opt::WeightOpt` doesn't check or clamp anything.
"""
function Weighted(array::AbsArray, vector::AbsVec = ones(lastlength(array)); norm=true)
    weights = clamp.(vector, 0, typemax(eltype(vector)))
    norm && rmul!(weights, 1/sum(weights))
    Weighted(array, weights, WeightOpt(norm=norm))
end

Weighted(array::AbsArray, lo::Real, hi::Real; norm=true) = Weighted(array, ones(lastlength(array)), lo, hi; norm=norm)

function Weighted(array::AbsArray, vector::AbsVec, lo::Real, hi::Real; norm=true)
    weights = clamp.(vector, 0, typemax(eltype(vector)))
    norm && rmul!(weights, 1/sum(weights))
    Weighted(clamp.(array,lo,hi), weights, WeightOpt(norm=norm, clamp=true, lo=lo, hi=hi))
end

Weighted() = Weighted(Float64[])

WeightedMatrix(x::AbsMat, rest... ; kw...) = Weighted(x, rest...; kw...)
WeightedMatrix(x::AbsVec, rest... ; kw...) = Weighted(reshape(x,:,1), rest...; kw...)
WeightedMatrix(x::AbsVec{T}) where {T} = Weighted(reshape(x,:,1), T[1], def_opt)

Base.vec(x::WeightedMatrix) = Weighted(vec(x.array), x.weights, x.opt)

"""
    hcat(x::WeightedMatrix, y, z...)
    x + y

These will `hcat` the arrays, and `vcat` the weights.
"""
function Base.hcat(x::Weighted, zz::Vararg{Union{Weighted, AbstractArray}})
    if length(x.array)==0
        if length(zz)==1 return zz[1]
        else return hcat(zz...)
        end
    end
    arr = hcat(x.array, [array(z) for z in zz]...)
    wei = vcat(x.weights, [sureweights(z) for z in zz]...)
    Weighted(arr, wei, x.opt)
end

sureweights(x::Weighted) = x.weights
sureweights(x) = begin n = lastlength(x); ones(n) ./ n end

Base.:+(x::Weighted, y::Union{AbsArray, Weighted}) = hcat(x,y)

"""
    λ * Π

For `Π::Weighted` this will scale `Π.weights` by `λ`.
Allows `0.7*Π1 + 0.3*Π2` since `+` means `hcat`.
Does not affect whether it is normalised.

    λ .* Π
    λ .+ Π

This instead scales, or shifts, `Π.array`,
and adjusts box constraints if applicable.
Arbitrary broadcasting is not supported, just these two!
"""
Base.:*(λ::Number, x::Weighted) = Weighted(x.array, λ .* x.weights, x.opt)
Base.:*(x::Weighted, λ::Number) = Weighted(x.array, λ .* x.weights, x.opt)

Base.broadcasted(::typeof(*), λ::Number, x::Weighted) = _scaleby(x, λ)
Base.broadcasted(::typeof(*), x::Weighted, λ::Number) = _scaleby(x, λ)

Base.broadcasted(::typeof(+), λ, x::Weighted) = _shiftby(x, λ)
Base.broadcasted(::typeof(+), x::Weighted, λ) = _shiftby(x, λ)
Base.broadcasted(::typeof(-), x::Weighted, λ) = _shiftby(x, -λ)

using Setfield

function _scaleby(x, λ)
    o = x.opt
    if x.opt.clamp
        o = set(o, @lens(_.lo), λ * o.lo)
        o = set(o, @lens(_.hi), λ * o.hi)
    end
    Weighted(λ .* x.array, x.weights, o)
end
function _shiftby(x, λ::Number)
    o = x.opt
    if x.opt.clamp
        o = set(o, @lens(_.lo), λ + o.lo)
        o = set(o, @lens(_.hi), λ + o.hi)
    end
    Weighted(λ .+ x.array, x.weights, o)
end
_shiftby(x, λ::AbstractVector) = Weighted(λ .+ x.array, x.weights, unclamp(x.opt))

wname(o::WeightOpt, s::Union{String,Symbol}) = set(o, @lens(_.wname), string(s))
aname(o::WeightOpt, s::Union{String,Symbol}) = set(o, @lens(_.aname), string(s))
addname(o::WeightOpt, s::Union{String,Symbol})  = set(o, @lens(_.aname), o.aname * string(s))
addlname(o::WeightOpt, s::Union{String,Symbol}) = set(o, @lens(_.aname), string(s) * o.aname)

##### Copying

Base.copy(x::Weighted) = Weighted(Base.copy(x.array), Base.copy(x.weights), x.opt)
Base.copyto!(x::Weighted, y::Weighted) = begin copyto!(x.array, y.array); copyto!(x.weights, y.weights); x end

wcopy!(x::Weighted, vec::AbsVec) = (x.weights .= clamp.(vec,0,Inf); x)
wcopy(x::Weighted, vec::AbsVec) = Weighted(x.array, clamp.(vec,0,Inf), x.opt) # only weights are copied!

##### Comparisons

Base.isapprox(x::Weighted, y::Weighted; kw...) =
    isapprox(array(x), array(y); kw...) && isapprox(weights(x), weights(y); kw...)

##### Extractors

array(x::Weighted) = x.array
array(x) = x
Base.parent(x::Weighted) = x.array

"""
    weights(x::Weighted)
If `x` is supposed to have normalised weights, this returns `normalise(x.weights)`
which also clamps them to be positive.
"""
# weights(x::Weighted) = x.opt.norm ? x.weights .* (1/sum(x.weights)) : x.weights ## sum(positive, x.weights) gives exotic problems
weights(x::Weighted) = x.opt.norm ? normalise(x.weights) : x.weights

positive(x) = max(x,zero(x))

totweight(x::Weighted) = sum(positive, x.weights)
maxweight(x::Weighted) = maximum(x.weights)

Base.getindex(x::Weighted, rr::AbsVec{Int}, c::Colon) = Weighted(x.array[rr,c], x.weights,  x.opt) ## get some dimensions
Base.getindex(x::Weighted, r::Colon, cc::AbsVec{Int}) = Weighted(x.array[r,cc], x.weights[cc], x.opt) ## get some points

using Lazy: @forward

@forward Weighted.array  Base.size, Base.ndims, Base.eltype, Base.getindex, Base.setindex!,
    Base.iterate, Base.lastindex,
    Base.minimum, Base.maximum, Base.extrema

lastlength(x) = size(x, ndims(x))

using StaticArrays

svecs(x::Weighted{<:Matrix{T}}) where T = reinterpret(SVector{size(x,1),T}, vec(x.array)) ## not type-stable
svecs(x::Weighted{<:Matrix{T}}, vald::Val{D}) where {T,D} = (@assert D==size(x,1); reinterpret(SVector{D,T}, vec(x.array)))

trace(x::Weighted) = x.opt.trace

##### Standardisers

clampdoc = """
    clamp!(x::Weighted)
    clamp(x::Weighted)
Always clamps weights to be positive, and if flag `clamp=true` is set in `x.opt`,
then clamps `x.array` using `lo,hi` from `x.opt`.

    clamp(x, lo, hi)
These alter the flag `x.opt.clamp` & then proceed.
"""
@doc clampdoc
function Base.clamp!(x::Weighted)
    clamp!(x.weights, 0, Inf)
    x.opt.clamp && clamp!(x.array, x.opt.lo, x.opt.hi)
    x
end
@doc clampdoc
function Base.clamp(x::Weighted)
    weights = clamp.(x.weights, 0, Inf)
    array = x.opt.clamp ? clamp.(x.array, x.opt.lo, x.opt.hi) : copy(x.array)
    Weighted(array, weights, x.opt)
end

# Base.clamp!(x::Weighted, lo::Real, hi::Real) = begin x.opt = clamp(x.opt, lo, hi); clamp!(x) end
Base.clamp(x::Weighted, lo::Real, hi::Real) =
    clamp!(Weighted(copy(x.array), copy(x.weights), Base.clamp(x.opt, lo, hi)))

function Base.clamp(o::WeightOpt, lo=0.0, hi=1.0)
    o = set(o, @lens(_.clamp), true)
    o = set(o, @lens(_.lo), lo)
    o = set(o, @lens(_.hi), hi)
    o
end
unclamp(o::WeightOpt) = set(o, @lens(_.clamp), false)

unclamp(x::Weighted) = Weighted(x.array, x.weights, unclamp(x.opt))

"""
    normalise(x) # with an s, NB!
    normalise!(x)
Ensures weights are positive and sum to 1.
On `x::Weighted`... the mutating form checks whether `norm=true` in `x.opt`;
the copying (!) form sets this flag first.

	unnormalise(x)
Just alters the flag `x.opt.norm`.
"""
function normalise(x::AbsVec{T}) where {T}
    itot = 1/sum(positive,x)
    clamp.(x, zero(T), typemax(T)) .* itot
end
@doc @doc(normalise)
normalise!(x::AbsVec) = begin clamp!(x,0,Inf); rmul!(x, 1/sum(x)) end

normalise!(x::Weighted) = begin x.opt.norm && normalise!(x.weights); x end

normalise(o::WeightOpt) = set(o, @lens(_.norm), true)
unnormalise(o::WeightOpt) = set(o, @lens(_.norm), false)

normalise(x::Weighted) = Weighted(copy(x.array), normalise(x.weights), normalise(x.opt))

unnormalise(x::Weighted) =  Weighted(copy(x.array), copy(x.weights), unnormalise(x.opt))

##### Flatten

"""
    flatten(x::Weighted)
Makes a vector `[ reshaped(x.array,:) ; x.weights ]`.

    flatcopy!(x, v::Vector)
Copies numbers from `v = flatten(x)` back into `x`.

    flatcopy(x, v)
Uses only `size(x)` to know what shape `::Weighted` to make out of `v` (and `x.opt` for details).
"""
flatten(x::Weighted) = vcat(vec(x.array), copy(x.weights))

@doc @doc(flatten)
function flatcopy!(x::Weighted, flat::AbsVec{<:Number})
    x.array[:] .= @view flat[1:length(x.array)]
    x.weights[:] .= @view flat[end-length(x.weights)+1:end]
    x
end

@doc @doc(flatten)
function flatcopy(x::Weighted, flat::AbsVec)
    array = reshape(flat[1:length(x.array)], size(x))
    weights = flat[end-length(x.weights)+1:end]
    Weighted(array, weights, x.opt)
end

using CatViews
"""
    CatView(x::Weighted)
Gives a vector-like view of both `x.array` and `x.weights`, the last `k` components are weights.
"""
CatViews.CatView(x::Weighted) = CatView(x.array, x.weights)

##### Pretty

function Base.show(io::IO, x::Weighted) ## compact, re-digestable
    print(io, "Weighted( ")
    show(io, array(x))

    print(io, ", ")
    show(io, weights(x))

    x.opt.clamp && print(io, ", ", x.opt.lo, ", ", x.opt.hi )
    x.opt.norm  || print(io, "; norm=false")
    print(io, " )")
end

function Base.show(io::IO, o::WeightOpt) ## compact, re-digestable
	print(io, "WeightOpt(")
	o.norm || print(io, "norm=false, ")
	o.clamp && print(io, "clamp=true, lo=",o.lo,", hi=",o.hi,", ")
    o.aname !== "θ" && print(io, "aname=",o.aname,", ")
    o.wname != "p(θ)" && print(io, "wname=",o.wname,", ")
    o.like && print(io, "like=true")
	print(io, ")")
end

function Base.show(io::IO, m::MIME"text/plain", x::Weighted) ## full
    ioc = IOContext(io, :compact => true, :limit => true)

    print(io, "Weighted ", summary(x.array))
    if x.opt.like
        println(io, ", likelihood ", x.opt.aname,":")
    elseif x.opt.clamp
        println(io, ", clamped ", round(x.opt.lo, digits=3), " ≦ ", x.opt.aname, " ≦ ", round(x.opt.hi, digits=3), ":")
    else
        println(io, ", of unclamped ", x.opt.aname,":")
    end
    if ndims(x.array)==1 && eltype(x.array) <: Number
        Base.print_array(ioc, x.array')  ## print vector as a row
    else
        Base.print_array(ioc, x.array)
    end

    length(x.weights)==0 && return nothing ## don't remember why
    if x.opt.norm
        println(io, "\nwith normalised weights ",x.opt.wname,", ", summary(x.weights), ":")
    else
        println(io, "\nwith weights ",x.opt.wname,", ", summary(x.weights), ":")
    end
    Base.print_array(ioc, weights(x)')

    if x.opt.norm && !isapprox(totweight(x), 1)
        print(io, "\nlazily normalised by weights(), raw total = ", round(totweight(x); digits=3), ".")
    end
end

function Base.summary(io::IO, x::Weighted)
    print(io, "Weighted ", summary(x.array))
    x.opt.clamp && print(io, " clamped to ", round(x.opt.lo, digits=3), " ... ", round(x.opt.hi, digits=3), ",")
    if x.opt.norm
        print(io, " with normalised weights ", summary(x.weights))
    else
        print(io, " with weights ", summary(x.weights))
    end
end
