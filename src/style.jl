include("colorscheme.jl")
using PyPlot, DrWatson
using3D()
export rdspl, axis_zoomin!, add_identifiers!, nice_arrow!, figx, figy

DrWatson._wsave(s, fig::Figure) = fig.savefig(s, dpi = 600, transparent = false)

"""
    rdspl(x, n = 3)
Round `x` with `n` sigdigits for display purposes.
"""
rdspl(x::Real, n = 3) = round(x, sigdigits=n)
rdspl(x::AbstractVector, n = 3) = Tuple((round.(Float64.(x); sigdigits=n)))

PyPlot.rc("lines", lw = 2.6)
PyPlot.rc("errorbar", capsize = 6)
PyPlot.rc("axes", grid = true)
PyPlot.rc("grid", color = "0.75", alpha = 0.75)
PyPlot.rc("axes", axisbelow=true) # makes grid behind other plotting elements

MATHFONTSIZE = 34
PyPlot.rc("font", size = 28) # set default fontsize
PyPlot.rc("xtick", labelsize = 24)
PyPlot.rc("ytick", labelsize = 24)
PyPlot.rc("axes", labelsize = MATHFONTSIZE)
PyPlot.rc("legend", fontsize = 30)
# PyPlot.rc("font", family = "Times New Roman") # Serif main font
PyPlot.rc("font", family = "DejaVu Sans") # sans main font
# PyPlot.rc("mathtext", rm = "sanserif", fontset="dejavusans") # sans math font
PyPlot.rc("mathtext", rm = "serif", fontset="cm") # serif math font

for z in ("x", "y")
    PyPlot.rc("$(z)tick.major", size = 7, width = 1.2)
    PyPlot.rc("$(z)tick.minor", size = 3, visible = false)
end

figx = 16 # default width correspoding to full text width
figy = 5  # default height corresponding to 1 row of plots
PyPlot.rc("figure", figsize = (figx, figy))
PyPlot.rc("savefig", dpi = 600, transparent = true, format = "png")

# set default color cycle
PyPlot.rc("axes", prop_cycle = matplotlib.cycler(color=COLORS))

if false # test font
    figure(figsize=(10,8)); plot(rand(10), label = "\$a=\\int_0^\\infty xdx\$")
    xlabel("x label"); ylabel("\$H_2\$"); legend(); tight_layout()
end

if false # test color scheme
    fig = figure(figsize = (20, 15)) # show colors
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(223)
    ax5 = subplot(224)
    lw = 60
    L = length(COLORSCHEME)
    for (i, c) in enumerate(COLORS)
        chsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
        ax1.plot([0, 1], [0, 0] .+ i, color = c, lw = lw)
        ax1.set_title("color")
        ax2.plot([0, 1], [0, 0] .+ i, color = string(chsv[3]), lw = lw)
        ax2.set_title("brightness")
        ax3.plot([0, 1], [0, 0] .+ i, color = string(chsv[2]), lw = lw)
        ax3.set_title("saturation")
        x = 0:0.05:5π
        ax4.plot(x, rand(1:3) .* cos.(x .+ i/2) .+ rand(length(x))/5; color=c, lw = 2)
        ax5.bar(collect(1:4) .+ (i-1)/L, 0.5rand(4) .+ 0.5, 1/L; color=c)
    end
    fig = tight_layout()
end

bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)

"`add_identifiers!(fig = gcf(), axs = fig.get_axes(); xloc, yloc)`"
function add_identifiers!(fig = gcf(), axs = fig.get_axes(); 
        xloc = 0.975, yloc = 1
    )
    bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)
    for (i, ax) in enumerate(axs)
        l = collect('a':'z')[i]
        if ax.name ≠ "3d"
            ax.text(xloc, yloc, "$(l)"; transform = ax.transAxes,
            bbox = bbox, zorder = 99, va = "top")
        end
    end
end

"""
    nice_arrow!(ax, xc, yc, xspan, yspan;
        style = "<->", tex = "", xo = 0.2, yo = -0.2, color = "C0"
    )
"""
function nice_arrow!(ax, xc, yc, xspan, yspan;
    style = "<->", tex = "", xo = 0.2, yo = -0.2, color = "C0")

    ax.annotate("";
        xy=(xc-xspan/2, yc - yspan/2), xycoords="data",
        xytext=(xc+xspan/2, yc + yspan/2), textcoords="data",
        arrowprops = Dict(
            :arrowstyle=>style,
            :connectionstyle=>"arc3",
            :lw=>1.5, :color => color),
        zorder = 99, color
    )
    if tex != ""
        ax.text(xc + xo, yc + yo, tex; va = :center, color)
    end
end

"""
    axis_zoomin!(zoomin, origin, zbox, rbox, co = "C0"; kwargs...)
Create a zoomin box connecting two axes, the `zoomin` and `origin`.
The zoom-in box is in the `origin` axis, while the `zoomin` axis is the
box. `rbox` is the enclosing box of the `zoomin` axes, while `zbox`
is the small "zoom-in" box of the `origin` axis. They must be in the form
((x1, y1), (x2, y2)). `co` is color `α` the alpha setting.
"""
function axis_zoomin!(zoomin, origin, zbox, rbox, co = "C0";
    connect_lines = true, lw = 2.0, α = 1.0, dir = :right, set_lims::Bool=true)
    # plot box in zoomin axis
    line, = zoomin.plot(
        [rbox[1][1], rbox[2][1], rbox[2][1], rbox[1][1], rbox[1][1]],
        [rbox[1][2], rbox[1][2], rbox[2][2], rbox[2][2], rbox[1][2]],
        color=co, lw = lw, alpha = α
    )
    line.set_clip_on(false)
    # plot box in origin axis
    line, = origin.plot(
        [zbox[1][1], zbox[2][1], zbox[2][1], zbox[1][1], zbox[1][1]],
        [zbox[1][2], zbox[1][2], zbox[2][2], zbox[2][2], zbox[1][2]],
        color=co, lw = lw, alpha = α
    )
    line.set_clip_on(false)

    # xyA is zoomin axis, xyB is origin axis
    if connect_lines
    if dir == :right
        for e in 1:2
            con = matplotlib.patches.ConnectionPatch(
            xyA = (rbox[1][1], rbox[e][2]), xyB=(zbox[2][1], zbox[e][2]),
            coordsA="data", coordsB="data",
            axesA = zoomin, axesB=origin, color=co, lw = lw, alpha = α)
            con.set_clip_on(false)
            zoomin.add_artist(con)
        end
    elseif dir == :top
        for e in 1:2
            con = matplotlib.patches.ConnectionPatch(
            xyA = (rbox[e][1], rbox[1][2]), xyB=(zbox[e][1], zbox[2][2]),
            coordsA="data", coordsB="data",
            axesA = zoomin, axesB=origin, color=co, lw = lw, alpha = α)
            con.set_clip_on(false)
            zoomin.add_artist(con)
        end
    end
    end
    # Set limits for zoom axis
    if set_lims
        zoomin.set_xlim(zbox[1][1], zbox[2][1])
        zoomin.set_ylim(zbox[1][2], zbox[2][2])
    end
    # remove axis for zoomin
    zoomin.axis("off")
end
