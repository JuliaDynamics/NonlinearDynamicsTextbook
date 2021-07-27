using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using DynamicalSystems, InteractiveDynamics
import GLMakie
using AbstractPlotting
using LinearAlgebra
using Statistics

# Code for 3D animation
ds = Systems.lorenz()
u0 = trajectory(ds, 10; Ttr = 100)[end]
dt = 0.01
T = 5.0
framerate = 30
systemtitle = "Lorenz63"

ds = Systems.towel()
u0 = trajectory(ds, 10; Ttr = 100)[end]
dt = 1
T = 15 # total evolove time
framerate = 2
systemtitle = "towelmap"

# %%
ε = 1e-6 # sphere radius
N = 4000 # initial conditions nmaking hypersphere

function uniform_sphere(N, u0, ε)
    V = typeof(u0)
    us = Vector{V}(undef, N)
    i = 1
    while i ≤ N
        x = rand(V) .- 0.5
        if norm(x) ≤ 1/(2*√3)
            us[i] = u0 .+ (ε .* x)
            i += 1
        end
    end
    return us
end

us = uniform_sphere(N, u0, ε)
# I must renormalize point coordinates to a sphere of radius 1 always,
# because of the internal convertion to 32-bit accuracy, small numbers are bad!
usnormed = [Point3f0(2((u0 .- u)./ε)...,) for u in us]
colors = [RGBAf0(abs.(u)..., 1) for u in usnormed]
plotus = Observable(usnormed)

# fig = Figure(resolution = (1200, 600))
fig = Figure(resolution = (800, 600))
ax3D = fig[1,1] = Axis3(fig)
scatter!(ax3D, plotus; markersize = 2000, color = colors, strokewidth = 0.5)
# ax2 = fig[1, 2] = Axis(fig; width = 400) # this axis has values of principal lengths
# ax2.title = "logarithm of principal axes"
ax3D.xticklabelsvisible = true
ax3D.yticklabelsvisible = true
ax3D.zticklabelsvisible = true
ax3D.azimuth = 4.44
ax3D.elevation = 0.38
ax3D.title = rpad(systemtitle *", t = $(0)", 40, ' ')

# function principals(v)
#     M = transpose(reshape(reinterpret(Float64, v), 3, :))
#     S = SVD(M).S # rotation matrix
#     return S
# end
# s = Observable(principals(usnormed))
# barplot!(ax2, s, color = to_color.(COLORSCHEME[1:3]))

display(fig)

# %%
ut = trajectory(ds, T, u0; dt)
evolution = [trajectory(ds, T, u; dt) for u in us]
t = 0:dt:T

# TODO: Use absolute, realsize of pointcloud

record(
        fig, joinpath(@__DIR__, "volume_growth_$(systemtitle).mp4"), 1:length(evolution[1]);
        framerate
        ) do i
    newus = [Point3f0(2((ut[i] .- e[i])./ε)...,) for e in evolution]
    plotus[] = newus
    # s[] = principals(newus)
    # autolimits!(ax2)
    # must set new limits centered on u
    maxε = maximum(norm(u) for u in newus)
    ax3D.limits = map(j -> (-maxε, +maxε), (1,2,3))
    ax3D.title = systemtitle *", t = $(t[i])"
    ax3D.elevation = 0.38
    ax3D.azimuth = 4.44
end
