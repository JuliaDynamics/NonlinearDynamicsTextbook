# %% show AC, MI for choosing optimal τ
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

ττ = 0:90
ac = autocor(w, ττ)
smi = selfmutualinfo(w, ττ)
smi ./= maximum(smi)

fig = figure(figsize = (figx/2, figy))
ax = gca()
ax.plot(ττ, ac; label = "AC")
ax.plot(ττ, smi, label = "SMI")
τ = estimate_delay(w, "mi_min")
ax.scatter(τ, smi[τ]; color = "C2", s = 100, zorder = 99)
ax.set_ylim(0,1.05)
ax.set_xticks(0:30:90)
ax.set_xlabel("\$\\tau\$", labelpad = -20)
ax.legend()
fig.tight_layout(;pad = 0.25)
wsave(plotsdir("mutualinfo_ac_tau"), fig)