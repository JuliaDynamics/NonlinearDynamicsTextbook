using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using InteractiveDynamics
using DynamicalSystems
using GLMakie

using OrdinaryDiffEq: Vern9
ds = Systems.henonheiles()

u0s = (
    [0.0, -0.25, 0.42, 0.0], # chaotic
    [0.0, 0.1, 0.5, 0.0], # quasiperiodic
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0], # periodic
)

diffeq = (alg = Vern9(), dtmax = 0.05)
idxs = (1, 2, 4)

fig, obs = interactive_evolution(
    ds, u0s; idxs, tail = 2000, diffeq, colors = COLORS[1:3],
)

ax = content(fig[1,1])
ax.xlabel = "q₁"
ax.ylabel = "q₂"
ax.zlabel = "v₂"
# main.figure[AbstractPlotting.Axis][:names, :axisnames] = ("q₁", "q₂", "v₂")


# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 15)
