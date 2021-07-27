# %% Roessler trajectory -> PSOS -> Lorenz map
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using DynamicalSystems, GLMakie

N = 10000.0

ds = Systems.roessler()
tr = trajectory(ds, N; Ttr = 100.0)

fig = Figure(resolution = (2100, 650), )
display(fig)

# Plot trajectory
trplot = fig[1,1] = Axis3(fig, scenekw = (camera = cam3d!, raw = false))
lines!(trplot, columns(tr)...; color = COLORS[5], linewidth = 2.0)
trplot.elevation = 0.4826990816987242
trplot.azimuth = 4.045530633326985
trplot.zticklabelsvisible=false
trplot.xticklabelsvisible=false
trplot.yticklabelsvisible=false

# Plot plane and section
o = Point3f0(-10, 0, 0)
w = Point3f0(25, 0, 25)
p = FRect3D(o, w)
a = RGBAf0(0,0,0,0)
c = to_color(COLORS[2])
img = GLMakie.Makie.ImagePattern([c a; a c]); # This throws an error if it shows, you can ignore that
mesh!(trplot, p; color = img)

psos = poincaresos(ds, (2, 0.0), N; Ttr = 100.0)
scatter!(trplot, columns(psos)...; color = COLORS[1], markersize = 4000)

# Plot section separately
psplot = fig[1, 2] = Axis(fig)
scatter!(psplot, psos[:, 1], psos[:, 3]; color = COLORS[1])
psplot.xlabel = "xₙ"
psplot.ylabel = "zₙ"

LS = 50
TS = 40

psplot.xticklabelsize = TS
psplot.yticklabelsize = TS
psplot.xlabelsize = LS
psplot.ylabelsize = LS

# Plot lorenz map
loplot = fig[1, 3] = Axis(fig)
scatter!(loplot, psos[1:end-1, 3], psos[2:end, 3]; color = COLORS[1])
loplot.xlabel = "zₙ"
loplot.ylabel = "zₙ₊₁"
loplot.xticklabelsize = TS
loplot.yticklabelsize = TS
loplot.xlabelsize = LS
loplot.ylabelsize = LS

trplot.zlabelsize = TS
trplot.xlabelsize = TS
trplot.ylabelsize = TS

display(fig)

# GLMakie.save(plotsdir("roessler_map.png"), fig)