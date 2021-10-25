# %% permutation entropy
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

using Random, Combinatorics, Statistics
Random.seed!(356)
o = 3
N = 8
x = rand(N)

fig = figure(figsize = (figx, figx/3.5))
ax1 = subplot2grid((1, 3), (0, 0))
ax1.plot(1:N, x, marker = "o", color = "C0", mfc = "C2", ms = 15, zorder = 99)
ax1.set_title("timeseries", size = 28)
ax1.spines["bottom"].set_position("center")
ax1.spines["right"].set_color("none")
ax1.spines["top"].set_color("none")
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel("\$t\$")
ax1.set_ylabel("\$x\$", rotation = 0)
ax1.xaxis.set_label_coords(1.0, 0.5)
ax1.yaxis.set_label_coords(-0.05, 0.9)
ax1.set_ylim(-0.1, 1.1)
# ax1.axis("off")

# add pattern indications
for (s, n) in zip(("2", "3", "3"), (2, 4, 6))
    xs = [n, n, n+o-1, n+o-1]
    ym = minimum(x[xs[1]:xs[end]]) - 0.1
    yd = ym - 0.02
    ys = [ym, yd, yd, ym] #.+ 0.05
    ax1.plot(xs, ys, color = "k", lw = 1)
    ax1.text(mean(xs), minimum(ys)-0.02, "#=$s", ha="center", va = "top")
end

# ax1.annotate("4", (2, -0.2), xytext = (4, -0.2), ha="center")

ax2 = subplot2grid((1, 3), (0, 1), colspan = 2)
p = permutations([0, 0.5, 1], o) |> collect
counts = [0, 1, 3, 1, 0, 1]
for (i, a) in enumerate(p)
    ax2.plot((1:o) .+ i*o .+ 1, a, marker = "o",
    color = "C0", mfc = "C2", ms = 10)
    ax2.text(o÷2 + i*o + 1, 1.2, "$i")
    ax2.text(o÷2 + i*o + 1, -0.5, "$(counts[i])")
end
ax2.text(o+1, 1.2, "#", ha = "right")
ax2.text(o+1, 0.5, "pattern", ha = "right")
ax2.text(o+1, -0.5, "count", ha = "right")
ax2.set_ylim(-1, 1.5)
ax2.set_xlim(0, ax2.get_xlim()[2])
ax2.set_title("relative amplitude permutations, \$d=$o\$", size = 28)
ax2.axis("off")

fig.tight_layout()
fig.subplots_adjust(top = 0.9, bottom = 0.02, right = 0.95, left = 0.05, wspace=0.1, hspace = 0.1)

wsave(plotsdir("6", "permentropy"), fig)
