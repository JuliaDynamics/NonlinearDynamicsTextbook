using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("theme.jl"))

using DynamicalSystems
using GLMakie

ds = Systems.lorenz()

u0s =  [[10,20,40.0] .+ rand(3) for _ in 1:6]

fig, obs = interactive_trajectory(
    ds, u0s; tail = 1000, colors = to_color.(COLORS)
)

ax3D = content(fig[1,1][1,1])
ax3D.elevation = 0.24269908169872415
ax3D.azimuth = 5.645530633326986

display(fig)

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
