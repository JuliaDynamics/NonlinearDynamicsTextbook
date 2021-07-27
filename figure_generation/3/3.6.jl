using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

a = 1.4; b = 0.3
he = Systems.henon(; a, b)
fig = figure()

x0, y0 = columns(trajectory(he, 5000, Ttr = 100))
x1 = x0
y1 = @. 1 - a*x0^2 + y0

x2 = b .* x1
y2 = y1
y3 = x2
x3 = y2

ax1 = subplot(141)
ax1.scatter(x0, y0, color = COLORS[1], s = 1, label = "initial")
ax1.text(0.1, 0.05, "initial", transform=ax1.transAxes)

αprev = 0.75
ax2 = subplot(142)
ax2.scatter(x0, y0, color = COLORS[1], s = 0.25, alpha = αprev)
ax2.scatter(x1, y1, color = COLORS[2], s = 1, label = "stretch/bend")
ax2.text(0.1, 0.05, "bend/stretch", transform=ax2.transAxes)
ax2.arrow(0,-0.19, 0, 0.8, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[5], length_includes_head=true)
ax2.arrow(-1.5, 0.0,  0, -0.5, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[5], length_includes_head=true)
ax2.arrow(1.5, 0.0,  0, -0.5, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[5], length_includes_head=true)

ax3 = subplot(143)
ax3.scatter(x1, y1, color = COLORS[2], s = 0.25, alpha = αprev)
ax3.scatter(x2, y2, color = COLORS[3], s = 1)
ax3.text(0.1, 0.05, "squeeze", transform=ax3.transAxes)
ax3.arrow(-1.5, 0.0,  1, 0, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[5], length_includes_head=true)
ax3.arrow(1.5, 0.0,  -1, 0, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[5], length_includes_head=true)


ax4 = subplot(144)
ax4.scatter(x2, y2, color = COLORS[3], s = 0.25, alpha = αprev)
ax4.scatter(x3, y3, color = COLORS[4], s = 1, label = "reflect")
ax4.text(0.1, 0.05, "reflect", transform=ax4.transAxes)

ax4.plot([-2,2], [-2,2], ls = "dashed", lw = 2, color = COLORS[5])

# ticks and stuff
ax1.set_ylabel("\$y\$")
for ax in (ax1, ax2, ax3, ax4)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks(-1.5:0.5:1.5)
    ax.set_xticks(-1.5:1:1.5)
    ax.grid("on")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_xlabel("\$x\$", labelpad = -10)
    # ax.legend(markerscale = 4, scatterpoints=1)
end
tight_layout(pad = 0.25)
subplots_adjust(left = 0.04, bottom = 0.1, wspace = 0.1, top = 0.9)
add_identifiers!(fig; xloc = 0)
wsave(plotsdir("stretchinghenon"), fig)
