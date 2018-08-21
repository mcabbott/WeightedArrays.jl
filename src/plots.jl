
export pplot, pplot!

PLOTSIZE = 400
ALPHA = 0.4

pointsize(x, sz=1) = 44 .* sqrt.(sz .* weights(x)) .+ 1

using RecipesBase

@recipe function f(x::Weighted{<:Matrix, <:Vector}; sz=1) ## exlude both 1D and Tracked
    if size(x,1)==1
        out = Weighted(vec(x.array), x.weights, x.opt) ## just to use 1D type signature
    elseif size(x,1)>=4
        out = wPCA(x,2)(x)
    else

        legend --> :false
        lab    --> " "
        if x.opt.clamp
            extr = [x.opt.lo - 0.05, x.opt.hi + 0.05]
            xlim --> extr
            ylim --> extr
            # zlim --> [x.opt.lo, x.opt.hi] ## leaves too much space, in GR
        end

        if endswith(x.opt.aname,"PC")
            size   --> (1.5*PLOTSIZE, PLOTSIZE)
            ylim   --> pcaylim(x, 1.5)
        elseif size(x,1)==2
            size   --> (1.1*PLOTSIZE, PLOTSIZE)
        elseif size(x,1)==3
            size   --> (1.3*PLOTSIZE, 1.2*PLOTSIZE)
        end

        seriestype := :scatter
        markeralpha --> ALPHA ## nickname alpha=0.2 doesn't overwrite this
        #markerstrokewidth --> 0
        markersize := pointsize(x, sz)

        if size(x,1)==2
            xaxis --> grsafe(x.opt.aname)*"_1"
            yaxis --> grsafe(x.opt.aname)*"_2"

            out = array(x)[1,:], array(x)[2,:]

        elseif size(x,1)==3
            xaxis --> xaxis_grbug(x.opt.aname)
            yaxis --> yaxis_grbug(x.opt.aname)
            zaxis --> grsafe(x.opt.aname)*"_3"

            zcolor --> array(x)[3,:]

            out = array(x)[1,:], array(x)[2,:], array(x)[3,:]
        end
    end
    out
end

# using Requires
# @require Plots begin

    ## Hack to work around bug:  https://github.com/JuliaPlots/Plots.jl/issues/743
    function xaxis_grbug(s::Union{Symbol,String})
        # Plots.backend_name() != :gr && return string(s)*"_1"
        return "        "*grsafe(s)*"\\_1                                       "*grsafe(s)*"\\_2"
    end
    function yaxis_grbug(s::Union{Symbol,String})
        # Plots.backend_name() != :gr && return string(s)*"_2"
        return grsafe(s)*"\\_3"
    end

    ## GR doesn't like unicode, but likes latex style greek
    function grsafe(s::Union{Symbol,String})
        # Plots.backend_name() != :gr && return string(s)
        return reduce(*, grsafe(c) for c in string(s))
    end

# end

function grsafe(c::Char)
    c=='θ' && return "\\theta"
    c=='ϕ' && return "\\phi"
    c=='σ' && return "\\sigma"
    c=='π' && return "\\pi"
    c=='_' && return "\\_"
    c=='(' && return "\\("
    c==')' && return "\\)"
    string(c)
end

function pcaylim(x::Weighted, ratio=1.5)
    mat = array(x)

    ex1 = [extrema(mat[1,:])...]
    diff1 = ex1[2] - ex1[1]

    ex2 = [extrema(mat[2,:])...]
    diff2 = ex2[2] - ex2[1]

    if diff1 > ratio * diff2
        ex2 .*= ( diff1 / (ratio * diff2) ) ## expand the yilm
    end

    return ex2
end

## 1D with shadow
@recipe function f(x::Weighted{<:Vector}, shadow="yes"; sz=0.8) ## sz doesn't get passed to here :(
    size   --> (1.2*PLOTSIZE, 0.8*PLOTSIZE)
    legend --> :false
    if x.opt.clamp && x.opt.hi < Inf
        xlim --> [x.opt.lo - 0.05, x.opt.hi + 0.05]
    end
    xaxis --> grsafe(x.opt.aname) #"\\theta"
    yaxis --> grsafe(x.opt.wname)
    ylim  --> [0, 1.4*maximum(weights(x)) ]
    yticks --> [0, round(maximum(weights(x)),digits=2)]

    if shadow != "only"
        @series begin ## plot the points
            seriestype := :scatter
            markersize := pointsize(x, sz)
            markeralpha --> ALPHA
            # markerstrokewidth --> 0

            array(x), weights(x)
        end
    end
    if shadow != "no"
        @series begin ## plot the shadow
            seriestype := :line
            lab   := ""
            fill  := 0
            color := :black
            alpha := 0.5*ALPHA

            shadowxy(array(x), weights(x))
        end
    end
end

function shadowxy(x::Vector, y::Vector, smooth=(maximum(x)-minimum(x))/100, numpoints=321; fixmax=true)
    @assert length(x)==length(y)
    xlist = range(minimum(x)-10*smooth, stop=maximum(x)+10*smooth, length=numpoints)
    ylist = zeros(numpoints)
    invtst = 1/(2*smooth^2)
    for i=1:numpoints
        for j=1:length(x)
            ylist[i] += y[j] * exp( -(xlist[i]-x[j])^2 * invtst )
        end
    end
    scale = fixmax ? maximum(y) / maximum(ylist) : 1
    return xlist, scale .* ylist
end

@recipe function f(x::Weighted, fun::Function)
    unique!(fun(x))
end

@recipe function f(x::Weighted, λ::Number)
    rmul!(copy(x), λ)
end

# @recipe function f(m::MultiModel; sz=1)
#     legend --> :topleft
#     for i=1:length(m)
#         @series begin
#             lab --> grsafe(getname(m,i)) * " : " * string(round(totweight(m.priors[i]),3))
#             m[i]
#         end
#     end
# end

using Requires

function __init__()
# function init_plot_req()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin

        """
            pplot(x::Weighted)
            pplot!(x)
        The first calls `sPCA`, the second calls `rPCA` to plot on the same axes.

            pplot(x, f)
            pplot(x, λ)
        Ditto, but first applies function `f`, or scales by `λ`.
        """
        pplot(x::Weighted, f::Function=identity; kw...) = Plots.plot(sPCA(unique!(f(x))); kw...)
        pplot(x::Weighted, λ::Number; kw...) = Plots.plot(rmul!(sPCA(x),λ); kw...)

        pplot!(x::Weighted, f::Function=identity; kw...) = Plots.plot!(rPCA(unique!(f(x))); kw...)
        pplot!(x::Weighted, λ::Number; kw...) = Plots.plot!(rmul!(rPCA(x),λ); kw...)

    end
end
