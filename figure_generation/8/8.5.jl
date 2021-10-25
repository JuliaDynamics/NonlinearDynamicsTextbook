using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using InteractiveDynamics, Random
using DynamicalBilliards
import GLMakie
using GLMakie: to_color, RGBf0, RGBAf0
using DynamicalSystems

fig = GLMakie.Figure(resolution = (figx*100*0.8, figy*100), fontsize = 32)
axhe = GLMakie.Axis3(fig[1,1])
axbi = GLMakie.Axis3(fig[1,2])
display(fig)

# %% Natural measure of henon map
ε = 0.001
x, y, p = begin
    ds = Systems.henon()
    X = trajectory(ds, 10^7; Ttr = 100)
    p, b = binhist(X, ε)
    x = [a[1] for a in b]
    y = [a[2] for a in b]
    x, y, p ./ maximum(p)
end

using GLMakie: Point3f0, Vec3f0, Rect3D

GLMakie.meshscatter!(axhe, vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(2ε, 2ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = clamp.(p, 0, 0.25), alpha = 0.5
)
# GLMakie.ylims!(axhe, 1.5 .* extrema(y)...)
# GLMakie.xlims!(axhe, extrema(x)...)
GLMakie.zlims!(axhe, 0, 0.2)
# GLMakie.xticks!(fig; )
axhe.zlabel = "ρ"

axhe.azimuth = 7.175530633326986
axhe.elevation = 0.53
axhe.zticklabelsvisible = false
axhe.zlabeloffset = 20
axhe.xticklabelsize = 26
axhe.yticklabelsize = 26
axhe.title = "Hénon map"

# %% same but for mushroom billiard
w = 0.2
bd = billiard_mushroom(1, w, 1, 0)

pc = Particle(-0.01, 0.2, sqrt(3)/2)
pr = Particle(0.0, 1.2, 0.0)
pr2 = Particle(0.0, 1.9, 0.0)

colors = [to_color(COLORS[i]) for i in (1,2,4)]
cmaps = [[InteractiveDynamics.darken_color(c, 2), c, InteractiveDynamics.lighten_color(c, 2)] for c in colors]
particles = [pc, pr, pr2]

ε = 0.02

for i in 1:3
    x, y, p = begin
        bmap, = boundarymap(particles[i], bd, 10^7)
        X = Dataset(bmap)
        p, b = ChaosTools.binhist(X, ε)
        x = [a[1]/2 for a in b]
        y = [a[2] for a in b]
        
        x, y, p ./ maximum(p)
    end

    GLMakie.meshscatter!(axbi, vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(ε, ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = clamp.(p, 0.25, 1), colormap = cmaps[i])
end
GLMakie.zlims!(axbi, 0, 1)

axbi.zlabel = "ρ"
axbi.zticklabelsvisible = false
axbi.xlabel = "ξ"
axbi.ylabel = "sin(θ)"
axbi.azimuth = 10.225530633326988
axbi.elevation = 0.53
axbi.zlabeloffset = 20
axbi.xticklabelsize = 26
axbi.yticklabelsize = 26
axbi.title = "mushroom billiard"

# %% Save
GLMakie.save(plotsdir("8", "natural_densities.png"), fig)
