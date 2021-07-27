# %% CCM illustration
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

ds = Systems.lorenz()
tr = trajectory(ds, 100; Ttr = 100)
x, y, z = columns(tr)

τx = estimate_delay(x, "mi_min")
τy = estimate_delay(y, "mi_min")
τz = estimate_delay(z, "mi_min")

X = embed(x, 3, τx)
Y = embed(y, 3, τy)
Z = embed(z, 3, τz)

using3D()
fig = figure()
axts = fig.add_subplot(1, 3, 1)
axx = fig.add_subplot(1,3,2; projection="3d")
axy = fig.add_subplot(1,3,3; projection="3d")

axts.plot(x .- 10)
axts.plot(y .+ 20; color = "C2")
axts.set_xlim(0, 1000)
axts.set_xticklabels([])
axts.set_yticklabels([])

axx.plot3D(columns(X)...; lw = 1)
axy.plot3D(columns(Y)...; lw = 1, color = "C2")
for s in (:x, :y, :z)
    f = Symbol(:set_, s, :ticklabels)
    @eval axx.$(f)([])
    @eval axy.$(f)([])
    g = Symbol(:set_, s, :lim)
    @eval axx.$(g)(-15, 15)
    @eval axy.$(g)(-20, 20)
end

# Axis pretty-fication
axx.dist = 8
axy.dist = 8
axx.elev = 20
axy.elev = 20
axx.text3D(-10, -10, 15, "\$M_x\$", size = 40)
axy.text3D(-10, -10, 20, "\$M_y\$", size = 40)
axts.set_ylabel("timeseries")

axx.text(15, 12, 21, "b"; bbox = bbox, zorder = 99, va = "top")
axy.text(17, 17, 27.5, "c"; bbox = bbox, zorder = 99, va = "top")
# axy.text(45, 40, 54, "c"; bbox = bbox, zorder = 99, va = "top")

add_identifiers!(fig)
fig.tight_layout(pad=0.35)
fig.subplots_adjust(wspace = 0.2)
wsave(plotsdir("ccm"), fig)