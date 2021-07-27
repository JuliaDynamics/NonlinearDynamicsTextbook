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
dt = 0.005
N = 100
framerate = 30
systemtitle = "Lorenz63"
azimuth = 2.24
elevation = 0.38

ds = Systems.towel()
u0 = trajectory(ds, 10; Ttr = 100)[end]
dt = 1
N = 15 # total steps to take
framerate = 2
systemtitle = "towelmap"

elevation = 0.5
azimuth = 6.02

# %%
tinteg = tangent_integrator(ds, 3; u0)
Y = get_deviations(tinteg)

# Now I define functions that make an ellipsoid given the point vectors
function ellipsoid_transform(Y)
    svddecomp = svd(Y)
    d = svddecomp.S
    s = d ./ maximum(d)
    Q = svddecomp.U * svddecomp.Vt # rotation matrix
    θx = atan(Q[3,2],Q[3,3])
    θy = atan(-Q[3,1], sqrt(Q[3,2]^2+Q[3,3]^2))
    θz = atan(Q[2,1], Q[1,1])
    M = AbstractPlotting.rotationmatrix_x(θx) *
        AbstractPlotting.rotationmatrix_y(θy) *
        AbstractPlotting.rotationmatrix_z(θz) *
        AbstractPlotting.scalematrix(Vec3f0(s...))
    return M, s
end
function transform_mesh(msh, mat4x4)
    pos_trans = Point3f0.(Ref(mat4x4) .* to_ndim.(Point4f0, spheremesh.position, 0))
    AbstractPlotting.GeometryBasics.Mesh(pos_trans, AbstractPlotting.GeometryBasics.faces(msh))
end
function ellipsoid_from_pointcloud(v, msh = spheremesh)
    M, s = ellipsoid_transform(v)
    return transform_mesh(msh, M), s
end
spheremesh = AbstractPlotting.GeometryBasics.mesh(Sphere(Point3f0(0), 1))

ellipsoid, s = ellipsoid_from_pointcloud(Y, spheremesh)
colors = [RGBAf0((abs.(u)./norm(u))..., 1) for u in spheremesh.position]

ellipsobs = Observable(ellipsoid)

# Setup figure and plot everything

# fig = Figure(resolution = (1200, 600))
fig = Figure(resolution = (1200, 600)); display(fig)
# Plot trajectory
t = Observable(0.0)
tstr = lift(x -> systemtitle *", t = $(round(x; digits = 5))", t)
axtr = Axis3(fig[1,1]; title = tstr, titlealign = :left)
tr = trajectory(ds, 100N*dt, u0; dt)
dotobs = Observable([u0])
trobs = Observable([u0])
if ds isa ContinuousDynamicalSystem
    lines!(axtr, tr.data; linewidth = 0.2, color = RGBAf0(0,0,0,0.5), transparency = true )
    lines!(axtr, trobs; linewidth = 4.0, color = to_color(JULIADYNAMICS_COLORS[1]), overdraw = true)
else
    scatter!(axtr, tr.data; markersize = 1000, color = RGBAf0(0,0,0,0.5), strokewidth=0, transparency = true )
    scatter!(axtr, trobs; markersize = 2000, color = to_color(JULIADYNAMICS_COLORS[1]), overdraw = true)
end
scatter!(axtr, dotobs; markersize = 4000, marker=:diamond,
color = to_color(JULIADYNAMICS_COLORS[3]), overdraw = true)
axtr.azimuth = azimuth
axtr.elevation = elevation

# Plot ellipsoid
ax3D = Axis3(fig[1, 2]; aspect = :data)
mesh!(ellipsobs, color = colors)

ax3D.xticklabelsvisible = false
ax3D.yticklabelsvisible = false
ax3D.zticklabelsvisible = false
ax3D.xlabelvisible = false
ax3D.ylabelvisible = false
ax3D.zlabelvisible = false
ax3D.azimuth = azimuth
ax3D.elevation = elevation

ax2 = Axis(fig[1, 3]; width = 300) # this axis has values of principal lengths
ax2.title = "ellipsoid axes"
tightlimits!(ax2, Bottom())
sobs = Observable(s)
barplot!(ax2, 1:3, sobs; color = to_color.(COLORSCHEME[1:3]))

display(fig)

# %%
record(
        fig, joinpath(@__DIR__, "volume_growth_$(systemtitle).mp4"), 1:N;
        framerate
        ) do i
    DynamicalSystems.step!(tinteg, dt)
    Y = get_deviations(tinteg)
    ellipsobs[], sobs[] = ellipsoid_from_pointcloud(Y)
    # s[] = principals(newus)
    # autolimits!(ax2)
    # must set new limits centered on u
    # maxε = maximum(norm(u) for u in newus)
    # maxε = maximum(s)
    # ax3D.limits = map(j -> (-maxε, +maxε), (1,2,3))
    # Update trajectory and title and elevation
    u = get_state(tinteg)
    dotobs[] = [u]
    trobs[] = push!(trobs[], u)
    t[] = tinteg.t
end
