
PLOTSIZE = 400
ALPHA = 0.4

pointsize(x, sz=1) = 44 .* sqrt.(sz .* weights(x)) .+ 1

using RecipesBase

@recipe function f(x::Weighted{<:Matrix, <:Vector}; sz=1) # exlude both 1D and Tracked
    if size(x,1)==1
        out = Weighted(vec(x.array), x.weights, x.opt) # just to use 1D type signature
    elseif size(x,1)>=4
        out = wPCA(x,2)(x)
    else

        legend --> :false
        label  --> " "
        if x.opt.clamp
            extr = [x.opt.lo - 0.05, x.opt.hi + 0.05]
            xlims --> extr
            ylims --> extr
            # zlim --> [x.opt.lo, x.opt.hi] # leaves too much space, in GR
        end

        if endswith(x.opt.aname,"PC")
            size   --> (1.5*PLOTSIZE, PLOTSIZE)
            ylims   --> pcaylim(x, 1.5)
        elseif size(x,1)==2
            size   --> (1.1*PLOTSIZE, PLOTSIZE)
        elseif size(x,1)==3
            size   --> (1.3*PLOTSIZE, 1.2*PLOTSIZE)
        end

        seriestype := :scatter
        markeralpha --> ALPHA # nickname alpha=0.2 doesn't overwrite this
        #markerstrokewidth --> 0
        markersize := pointsize(x, sz)

        if size(x,1)==2
            xaxis --> grsafe(x.opt.aname)*"_1"
            yaxis --> grsafe(x.opt.aname)*"_2"

            out = array(x)[1,:], array(x)[2,:]

        elseif size(x,1)==3
            xaxis --> grsafe(x.opt.aname)*"_1"
            yaxis --> grsafe(x.opt.aname)*"_2"
            zaxis --> grsafe(x.opt.aname)*"_3"

            marker_z --> array(x)[3,:]

            out = array(x)[1,:], array(x)[2,:], array(x)[3,:]
        end
    end
    out
end

@recipe f(x::Weighted{<:Base.ReshapedArray}; sz=1) = copy(x) # these were otherwise not caught
@recipe f(x::Weighted{<:Base.ReinterpretArray}; sz=1) = copy(x)

# GR doesn't like unicode, but likes latex style greek
function grsafe(s::Union{Symbol,String})
    # Plots.backend_name() != :gr && return string(s)
    return reduce(*, grsafe(c) for c in string(s))
end


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
        ex2 .*= ( diff1 / (ratio * diff2) ) # expand the yilm
    end

    return ex2
end

## 1D with shadow
@recipe function f(x::Weighted{<:Vector}, shadow="yes"; sz=0.8) # sz doesn't get passed to here :(
    size   --> (1.2*PLOTSIZE, 0.8*PLOTSIZE)
    legend --> :false
    if x.opt.clamp && x.opt.hi < Inf
        xlims --> [x.opt.lo - 0.05, x.opt.hi + 0.05]
    end
    xaxis --> grsafe(x.opt.aname) #"\\theta"
    yaxis --> grsafe(x.opt.wname)
    ylims  --> [0, 1.4*maximum(weights(x)) ]
    yticks --> [0, round(maximum(weights(x)),digits=2)]

    if shadow != "only"
        @series begin # plot the points
            seriestype := :scatter
            markersize := pointsize(x, sz)
            markeralpha --> ALPHA
            # markerstrokewidth --> 0

            array(x), weights(x)
        end
    end
    if shadow != "no"
        @series begin # plot the shadow
            seriestype := :line
            label := ""
            fill  := 0
            seriescolor := :black
            seriesalpha := 0.5*ALPHA

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
    unique(fun(x))
end
@recipe function f(fun::Function, x::Weighted)
    unique(fun(x))
end

@recipe function f(x::Weighted, λ::Number)
    λ .* x
end
@recipe function f(λ::Number, x::Weighted)
    λ .* x
end
