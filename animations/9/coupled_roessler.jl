using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using InteractiveDynamics
using DynamicalSystems, OrdinaryDiffEq
import GLMakie

a = 0.2
b = 0.2
c = 5.7
amu = 0.02
ω1 = 1.1 + amu
ω2 = 1.1 - amu
k1 = 0.02
k2 = k1

y_init = [0.1, 0.2, 0., 0.11, 0.19, 0.1]

diffeq = (alg = Tsit5(), adaptive = false, dt = 0.02, reltol = 1e-9, abstol = 1e-9)

ds = Systems.coupled_roessler(y_init; ω1, ω2, a, b, c, k1, k2)

ps = Dict(
    1 => 1.0:0.01:1.5,
    2 => 1.0:0.01:1.5,
    5 => 0:0.01:6,
    6 => 0:0.001:0.2,
    7 => 0:0.001:0.2,
)
pnames = Dict(
    1 => "ωx",
    2 => "ωy",
    5 => "c",
    6 => "kx",
    7 => "ky",
)

u0s =  [ds.u0]

fig, obs = interactive_evolution_timeseries(
    ds, u0s, ps; tail = 1000, diffeq, idxs = [1,4], pnames,
    lims = ((-12,12), (-12,12))
)

ax = GLMakie.content(fig[1,1])
ax.xlabel = "x1"
ax.ylabel = "y1"

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 30, sleep_time = 2)