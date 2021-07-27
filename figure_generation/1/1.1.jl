# phase space plots
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

tscolor = COLORS[1]

fig = figure(figsize = (figx, figx/2))

ax = subplot(231)

ds = Systems.standardmap(k=1)
ic = [[0.1, 0.1], [2.5, 0.4], [1.88, 3.25]]
xs = range(0, stop = 2π, length = 7);
ys = range(0, stop = 2π, length = 3);
ys = copy(xs)
iters = 2000
for (i, c) in enumerate(ic)
    tr = trajectory(ds, iters, c)
    ax.scatter(columns(tr)..., s = 5)
end
# ax.set_xticks([0, 2π])
ax.set_xticklabels([])
ax.set_yticklabels([])
# ax.set_yticks([0, 2π])
# ax.set_yticks([0, 2π])
ax.set_yticklabels([])
ax.set_xlabel("\$\\theta\$", labelpad = -5)
ax.set_xlim(0, 2π)
ax.set_ylim(0, 2π)
ax.set_ylabel("\$v\$",labelpad=0)
ax.set_title("Standard map")

ax = subplot(234)
n = 100
tr = trajectory(ds, n, ic[1])
ax.plot(0:n, tr[:, 2], c = tscolor, marker = "o", lw = 1.0)
ax.set_yticks([0, 2π])
ax.set_yticklabels(["0", "2\$\\pi\$"])
ax.set_ylabel("\$v\$",labelpad=-15)
ax.set_xticks([0, n])
ax.set_xlabel("\$n\$", labelpad = -20)

ds = Systems.lorenz()
using3D()
axlo3d = subplot(232, projection = "3d")
axlo3d.set_title("Lorenz-63", pad = 24)

llw = 2.0
tr = trajectory(ds, 100; Ttr = 100, dt = 0.01)
axlo3d.plot3D(columns(tr)..., color = COLORS[1], lw = llw/2)
tr = trajectory(ds, 10; Ttr = 1000, dt = 0.01)
axlo3d.plot3D(columns(tr)..., color = COLORS[2], lw = 2llw)
tr = trajectory(ds, 10, [10,20,40.0]; Ttr = 10, dt = 0.01)
axlo3d.plot3D(columns(tr)..., color = COLORS[3], lw = llw)
axlo3d.set_xlabel("\$x\$", labelpad = -10)
axlo3d.set_ylabel("\$y\$", labelpad = -10)
axlo3d.set_zlabel("\$z\$", labelpad = -10)
axlo3d.set_xticklabels([])
axlo3d.set_yticklabels([])
axlo3d.set_zticklabels([])
axlo3d.dist = 8

# ax.set_xticks([-20,20])
# ax.set_yticks([-20,20])
# ax.set_zticks([10,40])

tn = 30
axlo2d = subplot(235)
tr = trajectory(ds, tn; Ttr = 40, dt = 0.01)
axlo2d.plot(0:0.01:tn, tr[:, 1], c = tscolor)

axlo2d.set_yticks([-15,15])
axlo2d.set_ylabel("\$x\$", labelpad = -40)
axlo2d.set_xticks([0, tn])
axlo2d.set_xlabel("\$t\$", labelpad = -20)


ds = Systems.henonheiles()
u0s = (
    [0.0, -0.25, 0.42, 0.0], # chaotic
    [0.0, 0.1, 0.5, 0.0], # quasiperiodic
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0], # periodic
)

axhe3d = subplot(233, projection = "3d")
axhe3d.clear()
axhe3d.set_title("Hénon–Heiles", pad = 24)
for (i, u0) in enumerate(u0s)
    tr = trajectory(ds, 50, u0; Ttr = 50)
    axhe3d.plot3D(tr[:, 1], tr[:, 2], tr[:, 3], c = COLORS[i], alpha = 1,
    lw = i == 3 ? 4 : 2)
end

axhe3d.set_xticklabels([])
axhe3d.set_yticklabels([])
axhe3d.set_zticklabels([])
axhe3d.set_xlabel("\$x\$"; labelpad = -10)
axhe3d.set_ylabel("\$y\$"; labelpad = -10)
axhe3d.set_zlabel("\$v_x\$",labelpad=-10)
axhe3d.dist = 8

axhe2d = subplot(236)
n = 150
tr = trajectory(ds, n, [0.0, -0.25, 0.48, 0.0]; Ttr = 100)
axhe2d.plot(0:0.01:n, tr[:, 2], c = tscolor)
axhe2d.set_xticks([0, n])
axhe2d.set_yticks([-0.5, 0.5])
axhe2d.set_ylabel("\$y\$", labelpad = -40)
axhe2d.set_xlabel("\$t\$", labelpad = -20)


fig.tight_layout(pad = 0.3)
fig.subplots_adjust(right = 0.95)
# wsave(plotsdir("trajectories"), fig)
