struct CyclicContainer <: AbstractVector{String}
    c::Vector{String}
    n::Int
end
CyclicContainer(c) = CyclicContainer(c, 0)

Base.length(c::CyclicContainer) = length(c.c)
Base.size(c::CyclicContainer) = size(c.c)
Base.getindex(c::CyclicContainer, i) = c.c[mod1(i, length(c.c))]
function Base.getindex(c::CyclicContainer)
    c.n += 1
    c[c.n]
end
Base.iterate(c::CyclicContainer, i = 1) = iterate(c.c, i)

COLORS = [
    "#6F4AC7",
    "#33CBD8",
    "#1B1B1B",
    "#E82727",
    "#535D7F",
    "#A6D210",
]

COLORSCHEME = COLORS

CCOLORS = CyclicContainer(COLORS)
LINESTYLES = CyclicContainer(["-", ":", "--", "-."])

export COLORSCHEME, COLORS, CCOLORS, LINESTYLES