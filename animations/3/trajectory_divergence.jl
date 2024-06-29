using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using InteractiveDynamics, DynamicalSystems, GLMakie

using OrdinaryDiffEq: Tsit5
ds = Systems.lorenz()

u0 = [10,10,10.0]
u0s =  [u0 .+ i*1e-3 for i in 1:3]

diffeq = (alg = Tsit5(), dtmax = 0.01)

figure, obs = interactive_trajectory_timeseries(
    ds, u0s; tail = 1000, diffeq, colors = COLORS[1:3],
    linekwargs = (linewidth = 2.0,)
)

content(figure[1,1]).azimuth = 5.375530633327001
content(figure[1,1]).elevation = 0.3926990816987245

framerate = 30
total_time = 15 # in seconds
framen = framerate*total_time

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
