using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("theme.jl"))

using DynamicalSystems
using GLMakie

# Standard map trajectories
ds = Systems.standardmap(; k = 1.0)
u0s = [[θ, p] for θ ∈ 0:2π for p ∈ 0:2π]
lims = ((0, 2π), (0, 2π))

figure, obs = interactive_trajectory(
    ds, u0s; tail = 1000, lims,
)

# extract the main axis to add custom labels:
main = content(figure[1, 1][1,1])
main.xlabel = "θ"
main.ylabel = "p"

display(figure)

# %%
record_interaction(figure, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
