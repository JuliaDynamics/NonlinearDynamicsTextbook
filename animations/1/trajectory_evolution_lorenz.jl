using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using InteractiveDynamics
using DynamicalSystems
using GLMakie

using OrdinaryDiffEq: Tsit5
ds = Systems.lorenz()

u0s =  [[10,20,40.0] .+ rand(3) for _ in 1:6]

diffeq = (alg = Tsit5(), dtmax = 0.01)

fig, obs = interactive_evolution(
    ds, u0s; tail = 1000, diffeq, colors = to_color.(COLORS)
)

ax3D = content(fig[1,1])
ax3D.elevation = 0.24269908169872415
ax3D.azimuth = 5.645530633326986

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
