export COLORSCHEME, COLORS, MARKERS, LINESTYLES
export figuretitle!, axesgrid, subplotgrid
export label_axes!, space_out_legend!, textbox!
export lighten, invert_luminance

########################################################################################
# Colorscheme
########################################################################################
# You can pick colorschemes by setting the `ENV` variable `COLORSCHEME`
COLORSCHEMES = Dict(
    "JuliaDynamics" => [
        "#7143E0",
        "#191E44",
        "#0A9A84",
        "#AF9327",
        "#791457",
        "#6C768C",
    ],
    "JuliaDynamicsLight" => [ # for usage with dark background
        "#855DE4",
        "#B7BEF1",
        "#15C1A5",
        "#DCC261",
        "#EC59BB",
        "#737A8C",
    ],
    "Petrol" => [
        "#006269",
        "#BD5DAA",
        "#171B37",
        "#86612A",
        "#691105",
        "#00A9B5",
    ],
    "CloudySky" => [
        "#0099CC",
        "#67757E",
        "#1B1D4B",
        "#D07B17",
        "#6F0D4D",
        "#0D9276",
    ],
    "Flames" => [
        "#84150F",
        "#D65A35",
        "#E2B830",
        "#36454F",
        "#B2BEB5",
        "#9C278C",
    ],
    "GreenMetal" => [
        "#35B15A",
        "#748790",
        "#125D1D",
        "#BBA222",
        "#2B33AD",
        "#2B2931",
    ],
)

# ENV["COLORSCHEME"]  = "JuliaDynamicsLight" # change this to test
TEST_NEW_THEME = false
COLORSCHEME = COLORSCHEMES[get(ENV, "COLORSCHEME", "JuliaDynamics")]
BGCOLOR = get(ENV, "BGCOLOR", :transparent)
AXISCOLOR = get(ENV, "AXISCOLOR", :black)

mutable struct CyclicContainer{V} <: AbstractVector{V}
    c::Vector{V}
    n::Int
end
CyclicContainer(c) = CyclicContainer(c, 0)

Base.size(c::CyclicContainer) = size(c.c)
Base.length(c::CyclicContainer) = length(c.c)
Base.iterate(c::CyclicContainer, state=1) = Base.iterate(c.c, state)
Base.getindex(c::CyclicContainer, i::Int) = c.c[mod1(i, length(c))]
Base.getindex(c::CyclicContainer, i::AbstractVector) = getindex.(Ref(c.c), i)

function Base.getindex(c::CyclicContainer)
    c.n += 1
    c[c.n]
end

COLORS = CyclicContainer(COLORSCHEME)

########################################################################################
# Set Makie theme
########################################################################################
# The rest require `Makie` accessible in global scope
MARKERS = CyclicContainer([:circle, :dtriangle, :rect, :star5, :xcross, :diamond])
# Linestyles implement a better dash-dot than the original default (too much whitespace)
# and a second dashed style with longer lines between dashes
LINESTYLES = CyclicContainer([:solid, :dash, :dot, Linestyle([0, 3, 4, 5, 6]), Linestyle([0, 5, 6])])

cycle = Cycle([:color, :marker], covary = true)
_FONTSIZE = 16
_LABELSIZE = 20

default_theme = Makie.Theme(
    # Main theme (colors, markers, etc.)
    backgroundcolor = BGCOLOR,
    palette = (
        color = COLORSCHEME,
        marker = MARKERS,
        linestyle = LINESTYLES,
        patchcolor = COLORSCHEME,
    ),
    linewidth = 3.0,
    # Sizes of figure and font
    Figure = (
        size = (1000, 600),
        figure_padding = 20,
    ),
    fontsize = _FONTSIZE,
    Axis = (
        xlabelsize = _LABELSIZE,
        ylabelsize = _LABELSIZE,
        titlesize = _LABELSIZE,
        backgroundcolor = BGCOLOR,
        xgridcolor = AXISCOLOR,
        ygridcolor = AXISCOLOR,
        xtickcolor = AXISCOLOR,
        ytickcolor = AXISCOLOR,
        bottomspinecolor = AXISCOLOR,
        topspinecolor = AXISCOLOR,
        leftspinecolor = AXISCOLOR,
        rightspinecolor = AXISCOLOR,
        xlabelcolor = AXISCOLOR,
        ylabelcolor = AXISCOLOR,
        yticklabelcolor = AXISCOLOR,
        xticklabelcolor = AXISCOLOR,
        titlecolor = AXISCOLOR,
    ),
    Legend = (
        patchsize = (40f0, 20),
    ),
    Colorbar = (
        gridcolor = AXISCOLOR,
        tickcolor = AXISCOLOR,
        bottomspinecolor = AXISCOLOR,
        topspinecolor = AXISCOLOR,
        leftspinecolor = AXISCOLOR,
        rightspinecolor = AXISCOLOR,
        labelcolor = AXISCOLOR,
        ticklabelcolor = AXISCOLOR,
        titlecolor = AXISCOLOR,
    ),
    # This command makes the cycle of color and marker
    # co-vary at the same time in plots that use markers
    ScatterLines = (cycle = cycle, markersize = 5),
    Scatter = (cycle = cycle, markersize = 15),
    Band = (cycle = :color,),
    Lines = (cycle = Cycle([:color, :linestyle], covary = true),),
    Label = (textsize = _LABELSIZE,)
)

set_theme!(default_theme)

# Testing style (colorscheme)
if TEST_NEW_THEME
    # See also this website to see how the colorscheme looks for colorblind
    # https://davidmathlogic.com/colorblind
    using Random
    fig = Figure(size = (900, 600)) # show colors
    ax6 = Axis(fig[2,3])
    ax5 = Axis(fig[2,2])
    ax4 = Axis(fig[2,1])
    ax1 = Axis(fig[1,1]; title = "color")
    ax2 = Axis(fig[1,2]; title = "brightness")
    ax3 = Axis(fig[1,3]; title = "saturation")
    linewidth = 60
    L = length(COLORS)
    function graycolor(s)
        x = round(Int, 100s)
        return "gray"*string(x)
    end
    barpos = Random.shuffle(1:4L)
    for (i, c) in enumerate(COLORS)
        chsv = convert(Makie.Colors.HSV, to_color(c))
        lines!(ax1, [0, 1], [0, 0] .+ i; color = c, linewidth)
        lines!(ax2, [0, 1], [0, 0] .+ i; color = graycolor(chsv.v), linewidth)
        lines!(ax3, [0, 1], [0, 0] .+ i; color = graycolor(chsv.s), linewidth)
        local x = 0:0.05:5π
        lines!(ax4, x, rand(1:3) .* cos.(x .+ i/2) .+ rand(length(x))/5; color=c, linewidth = 4)
        barplot!(ax5, barpos[collect(1:4) .+ (i-1)*4], 0.5rand(4) .+ 0.5; width = 1, gap=0,color=c)
        scatterlines!(ax6, rand(3), rand(3); linewidth = 4, markersize = 30, color=c)
    end
    display(fig)
end


########################################################################################
# Convenience functions
########################################################################################
MakieText = Union{Symbol, <: AbstractString}

"""
    figuretitle!(fig, title; kwargs...)

Add a title to a `Figure`, that looks the same as the title of an `Axis`
by using the same default font. `kwargs` are propagated to `Label`.
"""
function figuretitle!(fig, title;
        valign = :bottom, padding = (0, 0, 0, 0),
        font = "TeX Gyre Heros Bold", # same font as Axis titles
        kwargs...,
    )
    Label(fig[0, :], title;
        tellheight = true, tellwidth = false, fontsize = _LABELSIZE, valign, padding, font, kwargs...
    )
    return
end

"""
    axesgrid(m, n; kwargs...) -> fig, axs

Create a grid of `m` rows and `n` columns of axes in a new figure
and return the figure and the `Matrix` of axis.

## Keyword arguments

- `sharex/sharey = false`: make every row share the y axis and/or every column
  share the x axis. In this case, tick labels are hidden from the shared axes.
- `titles::Vector{String}`: if given, they are used as titles for the axes of the top row.
  Can also be a single `String`, in which case it is used for all axes.
- `xlabels::Vector{String}`: if given, they are used as x labels of the axes
  in the bottom row. Can also be a single `String`, in which case it is used for all axes.
- `ylabels::Vector{String}`: if given, they are used as y labels of the axes in the
  leftmost column. Can also be a single `String`, in which case it is used for all axes.
- `title::String`: if given, it is used as super-title for the entire figure
  using the `figuretitle!` function.
- `kwargs...`: all further keywords are propagated to `Figure`.
"""
function axesgrid(m, n;
        sharex = false, sharey = false, titles = nothing,
        xlabels = nothing, ylabels = nothing, title = nothing, kwargs...
    )
    fig = Makie.Figure(; kwargs...)
    axs = Matrix{Axis}(undef, m, n)
    for i in 1:m
        for j in 1:n
            axs[i,j] = Axis(fig[i,j])
        end
    end
    if sharex
        for j in 1:n
            Makie.linkxaxes!(axs[:,j]...)
            for i in 1:m-1; Makie.hidexdecorations!(axs[i,j]; grid = false); end
        end
    end
    if sharey
        for i in 1:m # iterate through rows
            Makie.linkyaxes!(axs[i,:]...)
            for j in 2:n; Makie.hideydecorations!(axs[i,j]; grid = false); end
        end
    end
    if !isnothing(titles)
        for j in 1:n
            axs[1, j].title = titles isa MakieText ? titles : titles[j]
        end
    end
    if !isnothing(xlabels)
        for j in 1:n
            axs[end, j].xlabel = xlabels isa MakieText ? xlabels : xlabels[j]
        end
    end
    if !isnothing(ylabels)
        for i in 1:m
            axs[i, 1].ylabel = ylabels isa MakieText ? ylabels : ylabels[i]
        end
    end
    !isnothing(title) && figuretitle!(fig, title)
    return fig, axs
end
const subplotgrid = axesgrid


"""
    label_axes!(axs::Array{Axis};
        valign = :top, halign = :right, pad = 5,
        labels = range('a'; step = 1, length = length(axs)),
        add_box = false, boxkw = NamedTuple(), kw...
    )

Add labels (like a,b,c,...) to all axes.
Keywords customly adjust location, and `kw` are propagated to `Label`.
If chosen, a box is added around the label with options `boxkw` propagated to `Box`.
"""
function label_axes!(axs;
        labels = range('a'; step = 1, length = length(axs)),
        transformation = x -> "("*string(x)*")",
        valign = :top, halign = :right,
        pad = 5, add_box = false, kwargs...,
    )

    lbs = @. string(transformation(labels))
    # Create padding from alignment options
    padding = [0,0,0,0]
    if halign == :right
        padding[2] = pad
    elseif halign == :left
        padding[1] = pad
    end
    if valign == :top
        padding[3] = pad
    elseif valign == :bottom
        padding[4] = pad
    end

    for (i, ax) in enumerate(axs)
        @assert ax isa Axis
        gc = ax.layoutobservables.gridcontent[]
        x = gc.parent[gc.span.rows, gc.span.cols]
        # Currently `Label` has no way of having a box around it
        lab = Label(x, lbs[i];
            tellwidth=false, tellheight=false,
            valign, halign, padding, font = :bold, justification = :center,
            kwargs...
        )
        if add_box
            # but we can access the internals and get the box of the label,
            # and then make an actual box around it
            bx = Box(first(axs).parent; bbox = lab.layoutobservables.computedbbox, color = "transparent")
            Makie.translate!(bx.blockscene, 0, 0, -1)
        end

        # TODO: This is a much better option:
        # gc = ax.layoutobservables.gridcontent[]
        # pos = gc.parent[gc.span.rows, gc.span.cols]
        # Textbox(pos;
        #     placeholder = string(round(rmi; sigdigits = 2)),
        #     textcolor_placeholder = :black, valign = :top, halign = :right,
        #     tellwidth = false, tellheight=false, boxcolor = (:white, 0.75),
        #     textpadding = (4, 4, 4, 4)
        # )

    end
    return
end

"""
    textbox!(ax::Axis, text::AbstractString; kw...)

Add a small textbox to `ax` containing the given `text`.
By default, the textbox is placed at the top right corner with proper alignment,
and it is slightly transparent. See the source code for the keywords you need
to adjust for different placement or styling.
"""
function textbox!(ax, text; kw...)
    gc = ax.layoutobservables.gridcontent[]
    pos = gc.parent[gc.span.rows, gc.span.cols]
    Textbox(pos;
        placeholder = text,
        textcolor_placeholder = :black, valign = :top, halign = :right,
        tellwidth = false, tellheight=false, boxcolor = (:white, 0.75),
        textpadding = (2, 2, 2, 2), kw...
    )
    return
end

"""
    space_out_legend!(legend)

Space out the contents of a given legend, so that the banks are spaced equidistantly
and cover the full width available for the legend. This function is supposed to be
called for horizontal legends that should span the full width of a column
and hence are placed either on top or below an axis.
"""
function space_out_legend!(legend)
    Makie.update_state_before_display!(legend.parent)
    # ensure that the width and height are told correctly
    legend.tellwidth[] = false
    legend.tellheight[] = true
    # update axis limits etc, the legend adjustment must come last
    w_available = legend.layoutobservables.suggestedbbox[].widths[1]
    w_used = legend.layoutobservables.computedbbox[].widths[1]
    difference = w_available - w_used
    legend.colgap[] += difference / (legend.nbanks[] - 1)
    return
end
########################################################################################
# Color manipulation
########################################################################################
"""
    invert_luminance(color)

Return a color with same hue and saturation but luminance inverted.
"""
function invert_luminance(color)
    c = to_color(color)
    hsl = Makie.HSLA(c)
    l = 1 - hsl.l
    neg = Makie.RGBA(Makie.HSLA(hsl.h, hsl.s, l, hsl.alpha))
    return neg
end

"""
    lighten(c, f = 1.2)

Lighten given color `c` by multiplying its luminance with `f`.
If `f` is less than 1, the color is darkened.
"""
function lighten(c, f = 1.2)
    c = to_color(c)
    hsl = Makie.HSLA(c)
    neg = Makie.RGBAf(Makie.HSLA(hsl.h, hsl.s, clamp(hsl.l*f, 0.0, 1.0), hsl.alpha))
    neg = Makie.RGBf(Makie.HSL(hsl.h, hsl.s, clamp(hsl.l*f, 0.0, 1.0)))
    return neg
end

########################################################################################
# I/O
########################################################################################
if isdefined(Main, :DrWatson)
    # Extension of DrWatson's save functionality for default CairoMakie saving
    function DrWatson._wsave(filename, fig::Makie.Figure, args...; kwargs...)
        if isdefined(Main, :CairoMakie)
            # Always save via CairoMakie for high quality if possible
            CairoMakie.save(filename, fig, args...; px_per_unit = 2, kwargs...)
        else
            Makie.save(filename, fig, args...; kwargs...)
        end
    end

    # Using FileIO's load to make figures for black slides
    """
        negate_remove_bg(file; threshold = 0.02, bg = :white, overwrite = false)

    Create an negated version of the image at `file` with background removed,
    so that it may be used in environments with dark background.
    The `threshold` decides when a pixel should be made transparent.
    If the image already has a dark background, pass `bg = :black` instead,
    which will not negate the image but still remove the background.
    """
    function negate_remove_bg(file; threshold = 0.02, bg = :white, overwrite = false)
        img = DrWatson.FileIO.load(file)
        x = map(px -> invert_color(px, bg, threshold), img)
        if overwrite
            newname = file
        else
            name, ext = splitext(file)
            newname = name*"_inv"*ext
        end
        DrWatson.FileIO.save(newname, x)
    end

    function invert_color(px, bg = :white, threshold = 0.02)
        hsl = Makie.HSLA(to_color(px))
        l = (bg == :white) ? (1 - hsl.l) : hsl.l
        neg = Makie.RGB(Makie.HSL(hsl.h, hsl.s, l))
        α = abs2(neg) < threshold ? 0 : hsl.alpha
        Makie.RGBA(neg, α)
    end
    function negate_remove_save(filename, fig::Makie.Figure)
        DrWatson.wsave(filename, fig)
        negate_remove_bg(filename; overwrite = true)
    end

    """
        remove_bg(file; threshold = 0.02, overwrite = false)

    Remove the background for figure in `file` (all pixels with luminosity > 1 - threshold).
    Either overwrite original file or make a copy with suffix _bgr.
    """
    function remove_bg(file; threshold = 0.02, overwrite = false)
        img = DrWatson.FileIO.load(file)
        x = map(px -> make_transparent_pixel(px, threshold), img)
        if overwrite
            newname = file
        else
            name, ext = splitext(file)
            newname = name*"_bgr"*ext
        end
        DrWatson.FileIO.save(newname, x)
    end
    function make_transparent_pixel(px, threshold = 0.02)
        α = Makie.RGBA(to_color(px)).alpha
        c = Makie.RGB(to_color(px))
        α = abs2(c) > 1 - threshold ? 0 : α
        return Makie.RGBA(c, α)
    end

end
