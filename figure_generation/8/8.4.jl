# %% Standard map island shrinking
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot
grid = 0.0:0.5:2π
sm = Systems.standardmap()
ks = [0.1, 1.0, 4.0]

fig, axs = subplots(1,3; sharey = true)

for i in 1:3
    set_parameter!(sm, 1, ks[i])
    ax = axs[i]
    for θ in grid, p in grid
        u0 = SVector(θ, p)
        tr = trajectory(sm, 1000, u0)
        c = DynamicalSystems.lyapunovspectrum(sm, 4000, 1; u0 = u0)[1] > 0.01 ? "C0" : "C2"
        ax.plot(columns(tr)...; ls = "None", marker = "o", ms = 0.2, c = c)
    end
    ax.set_xlim(0, 2π)
    ax.set_ylim(0, 2π)
    ax.text(4.5, 5.5, "\$k=$(ks[i])\$"; bbox = bbox)
    ax.set_xticks(0:2:6)
end

axs[1].set_ylabel("\$p\$", labelpad = -5)
axs[2].set_xlabel("\$\\theta\$", labelpad = -15)
fig.tight_layout(pad=0.3)
# fig.subplots_adjust(left = 0.05, bottom = 0.14, right = 0.95, top = 0.97, wspace = 0.06)
wsave(plotsdir("8","kam_sm"), fig)