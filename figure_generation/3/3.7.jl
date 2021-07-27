# %% chaoticity map
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

using OrdinaryDiffEq
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

fig = figure(figsize = (figx/2, figy))
ax = gca()
hh = Systems.henonheiles()
ics = Systems.henonheiles_ics(0.13, 15)
cmap = matplotlib.cm.get_cmap("viridis")

scres = nothing
for ic in ics
    psos = poincaresos(hh, (1, 0.0), 2000; u0 = ic)
    λ = lyapunov(hh, 10000; u0 = ic, diffeq...)
    λmax = 0.06
    v = clamp(λ/λmax, 0, 1)
    global scres = ax.scatter(psos[:, 2], psos[:, 4];
        color = cmap(v), s = 5, edgecolors = "0.5", linewidths = 0.1,
        vmin = 0, vmax = 1,
    )
end

cb = colorbar(scres)
cb.set_ticks([0, 1])
cb.set_ticklabels(["0", "0.06"])
cb.set_label("chaoticity (\$\\lambda_1\$)"; labelpad=-50)
ax.set_xticks(-0.5:0.4:0.8)
ax.set_xlabel("\$y\$", labelpad = -20)
ax.set_yticks(-0.5:0.3:0.5)
ax.set_ylabel("\$p_y\$", labelpad = -10)
fig.tight_layout(; pad = 0.25)
fig.subplots_adjust(bottom = 0.16, top = 0.95, left = 0.15, wspace = 0.01)
wsave(plotsdir("chaoticity"), fig)
