# uses file produced by `11.6_produce.jl`.
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using Random, OrdinaryDiffEq
import GLMakie

data = wload(datadir("FitzHugh", "spiralwave.jld2"))
@unpack uout, tvec, params = data

heatobs = GLMakie.Observable(uout[1][1,:,:])
tobs = GLMakie.Observable(0.0)
titobs = GLMakie.lift(t -> "t = $(t)", tobs)

fig = GLMakie.Figure(resolution = (600, 550))
ax = fig[1,1] = GLMakie.Axis(fig; title = titobs)
hmap = GLMakie.heatmap!(ax, heatobs; colormap = :tokyo, colorrange = (-0.2, 1))
cb = GLMakie.Colorbar(fig[1, 2], hmap; width = 20)
display(fig)

GLMakie.record(
        fig, projectdir("animations", "11", "fitzhugh_spiralwave.mp4"), 
        1:length(tvec); framerate = 15
    ) do i
    
    tobs[] = tvec[i]
    heatobs[] = uout[i][1,:,:]
end
