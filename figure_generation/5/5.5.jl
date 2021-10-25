# Henon gridding for dimension explanation
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

fig = figure()
ax = subplot(1,2,1)
ax2 = subplot(1,2,2)
he = Systems.henon()
N = 10000
tr = trajectory(he, N, Ttr = 100)
es = [0.2, 0.05, 0.01]
ns = zero(es)
ax.grid(false)
mini = minima(tr)

rectcols = ["C0", "C1", "C3"]

for (i, e) ∈ enumerate(es)
    if i ≤ 2
        for x in mini[1]:e:1.4
            ax.axvline(x; color = "k", lw = 0.5, alpha = 1/(3^i))
        end
        for y in mini[2]:e:0.4
            ax.axhline(y; color = "k", lw = 0.5, alpha = 1/(3^i))
        end
    end
    p, bins = binhist(tr, e)
    for b in bins
        r = matplotlib.patches.Rectangle(b, e, e; color = rectcols[i], ec = "k", lw = 1/i)
        ax.add_artist(r)
    end
    ns[i] = length(probabilities(tr, e))
end

ax.plot(tr[:, 1], tr[:, 2], ls = "None", marker = ".", color = COLORS[3], ms = 0.75, zorder = 99)
ax.set_yticks([])
ax.set_ylabel("\$y\$")
ax.set_xticks([])
ax.set_xlabel("\$x\$")

ax2.scatter(log.(1 ./ es), log.(ns); c = rectcols, s = 200, zorder = 99)
s = linreg(log.(1 ./ es), log.(ns))[2]
ax2.plot(log.(1 ./ es), log.(1 ./ es) .* s .+ 2, color = "C2")
ax2.text(3.5, 6, "\$\\Delta = $(rdspl(s))\$"; color = "C2")
ax2.set_xlabel("\$\\log ( 1/\\varepsilon)\$"; labelpad = -15)
ax2.set_ylabel("\$\\log ( M)\$")
ax2.set_xticks([2, 3, 4])
ax2.set_xticklabels([2, "", 4])
fig.tight_layout(pad = 0.3)
fig.subplots_adjust(wspace = 0.2)
wsave(plotsdir("5", "henon_gridding"), fig)
