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
    "#6D44D0",
    "#2CB3BF",
    "#1B1B1B",
    "#DA5210",
    "#866373",
    "#03502A"
]

COLORSCHEME = COLORS

CCOLORS = CyclicContainer(COLORS)
NONBLACK = CyclicContainer(COLORS[[1,2,4,5]])
BLACK = COLORS[3]

LINESTYLES = CyclicContainer(["-", "--", "-.",  ":",])

figx = 16 # default width correspoding to full text width, in inches
figy = 5  # default height corresponding to 1 row of plots, in inches

export COLORSCHEME, COLORS, CCOLORS, LINESTYLES, figx, figy