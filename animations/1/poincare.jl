using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("theme.jl"))

using DynamicalSystems
import GLMakie; GLMakie.activate!()

ds = Systems.henonheiles()

u0s = [
    [0.0, -0.25, 0.42, 0.0], # chaotic
    [0.0, 0.1, 0.5, 0.0], # quasiperiodic
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0], # periodic
]
trs = [trajectory(ds, 10_000, u0; Dt = 0.1)[1][:, SVector(1,2,3)] for u0 ∈ u0s]

fig = interactive_poincaresos_scan(trs, 2; linekw = (transparency = true,))

# %% use this to record it
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
