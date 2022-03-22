# Logistic map cobweb and timeseries
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveDynamics

lo = Systems.logistic(0.4; r=rrange[1])
fkwargs = [(linewidth = 4.0, color = COLORSCHEME[i+1]) for i in 1:3]

rrange = (rc = 1 + sqrt(8); [rc - 1e-3, rc - 1e-5, rc])
fig = interactive_cobweb(lo, rrange, 3; fkwargs)
record_interaction(
    fig, projectdir("animations", "4", "logistic_intermittency.mp4"); 
    total_time=15, sleep_time = 2
)
