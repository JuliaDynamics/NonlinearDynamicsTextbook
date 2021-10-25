# %% Roessler trajectory -> PSOS -> Lorenz map
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using DynamicalSystems
import GLMakie

N = 10000.0

ds = Systems.roessler()
tr = trajectory(ds, N; Ttr = 100.0)

fig = GLMakie.Figure(resolution = (2100, 650), )
display(fig)

# Plot trajectory
trplot = fig[1,1] = GLMakie.Axis3(fig, scenekw = (camera = GLMakie.cam3d!, raw = false))
GLMakie.lines!(trplot, columns(tr)...; color = COLORS[3], linewidth = 2.0)
trplot.elevation = 0.4826990816987242
trplot.azimuth = 4.045530633326985
trplot.zticklabelsvisible=false
trplot.xticklabelsvisible=false
trplot.yticklabelsvisible=false

# Plot plane and section
o = GLMakie.Point3f0(-10, 0, 0)
w = GLMakie.Point3f0(25, 0, 25)
p = GLMakie.FRect3D(o, w)
a = GLMakie.RGBAf0(0,0,0,0)
c = GLMakie.to_color(COLORS[4])
img = GLMakie.Makie.ImagePattern([c a; a c]); # This throws an error if it shows, you can ignore that
GLMakie.mesh!(trplot, p; color = img)

psos = poincaresos(ds, (2, 0.0), N; Ttr = 100.0)
GLMakie.scatter!(trplot, columns(psos)...; color = COLORS[1], markersize = 5000)

# Plot section separately
psplot = GLMakie.Axis(fig[1, 2])
GLMakie.scatter!(psplot, psos[:, 1], psos[:, 3]; color = COLORS[1])
psplot.xlabel = "xₙ"
psplot.ylabel = "zₙ"

LS = 50
TS = 40

psplot.xticklabelsize = TS
psplot.yticklabelsize = TS
psplot.xlabelsize = LS
psplot.ylabelsize = LS

# Plot lorenz map
loplot = fig[1, 3] = GLMakie.Axis(fig)
GLMakie.scatter!(loplot, psos[1:end-1, 3], psos[2:end, 3]; color = COLORS[1])
loplot.xlabel = "zₙ"
loplot.ylabel = "zₙ₊₁"
loplot.xticklabelsize = TS
loplot.yticklabelsize = TS
loplot.xlabelsize = LS
loplot.ylabelsize = LS

trplot.zlabelsize = TS
trplot.xlabelsize = TS
trplot.ylabelsize = TS
trplot.xlabeloffset = 20
trplot.ylabeloffset = 20
trplot.zlabeloffset = 20

# GLMakie.Label(fig[1, 1, GLMakie.TopRight()], "b");

display(fig)

GLMakie.save(plotsdir("4", "roessler_map.png"), fig)