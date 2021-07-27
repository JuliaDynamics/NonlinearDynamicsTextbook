using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using DynamicalSystems, InteractiveDynamics
using OrdinaryDiffEq
import GLMakie

ds = Systems.henonheiles()
diffeq = (alg = Vern9(),)
u0s = [
    [0.0, -0.25, 0.42, 0.0], # chaotic
    [0.0, 0.1, 0.5, 0.0], # quasiperiodic
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0], # periodic
]
trs = [trajectory(ds, 1000, u0; diffeq...)[:, SVector(1,2,3)] for u0 ∈ u0s]

fig, ax3D, ax2D = brainscan_poincaresos(trs, 2; linekw = (transparency = false,))

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
