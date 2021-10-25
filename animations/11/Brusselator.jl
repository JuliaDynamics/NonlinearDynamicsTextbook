# This file needs `src/brusselator_produce.jl` to have run with appropriate parameters.
# The `saveat` parameter of that file is the time output of the frames.
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using Random, OrdinaryDiffEq
import GLMakie

a = 9.0
b = 10.2
L = 50.0 # size of domain
N = 250  # no. of grid points
d = 1.7

prefix = "brusselator2D_dense"
fname = savename(prefix, @strdict(a, b, d, L, N, ic), "jld2")
data = load(datadir("Brusselator", fname))
uout = data["u"]
tvec = data["t"]

heatobs = GLMakie.Observable(uout[1][1,:,:])
tobs = GLMakie.Observable(0.0)
titobs = GLMakie.lift(t -> "t = $(t)", tobs)

fig = GLMakie.Figure(resolution = (600, 550))
ax = fig[1,1] = GLMakie.Axis(fig; title = titobs)
hmap = GLMakie.heatmap!(ax, heatobs; colormap = :inferno, colorrange = (0, 1.5))
cb = GLMakie.Colorbar(fig[1, 2], hmap; width = 20)
display(fig)

GLMakie.record(
        fig, projectdir("animations", "11", "brusselator_d=$(d)_b=$(b).mp4"), 
        1:length(tvec); framerate = 15
    ) do i
    
    tobs[] = tvec[i]
    heatobs[] = uout[i][1,:,:]
end
