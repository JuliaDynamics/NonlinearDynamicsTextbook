using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, InteractiveDynamics
import GLMakie
using AbstractPlotting.MakieLayout

# %% Sensitive dependence demonstration for Lorenz-63
using OrdinaryDiffEq: Tsit5, Vern9
ds = Systems.lorenz()
u0 = 3ones(dimension(ds))
u0s =  [u0 .+ 1e-3i for i in 1:3]

diffeq = (alg = Vern9(), dtmax = 0.01)

figure, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 100, diffeq, colors = COLORS,
)


# %% stretching and folding for logistic
lo = Systems.logistic()
X = trajectory(lo, 1000, Ttr = 100)
# close("all")
fig = figure()
# ax = subplot2grid((1,3), (0,0))
ax = subplot(121)
ax.scatter(X[1:end-1], X[2:end], color = COLORS[1], s = 5, label = "\$x_{i+2}\$")
ax.set_xlabel("\$x_{i}\$")
ax.set_yticks([0, 0.5, 1])
ax.scatter(X[1:end-2], X[3:end], color = COLORS[2], s = 5, label = "\$x_{i+2}\$")
ax.legend(markerscale = 5)

using3D()
# ax2 = subplot2grid((1,3), (0,1), projection = "3d", colspan = 2)
ax2 = subplot(122, projection = "3d")
ax2.scatter(X[1:end-2], X[2:end-1], X[3:end], color = COLORS[1], s = 5)

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_zticklabels([])
ax2.set_xlabel("\$x_{i}\$", labelpad = -10)
ax2.set_ylabel("\$x_{i+1}\$", labelpad = -10)
ax2.set_zlabel("\$x_{i+2}\$", labelpad = -10)

# fold arrow
s = (0.15, 0.5, 1.2)
e = (0.85, 0.5, 1.2)
ax2.quiver3D(e..., (s .- e)..., color = COLORS[3])
ax2.text3D(0.5, 0.1, 1.2, "Fold", color = COLORS[3])

# stretch arrow
s = (0.75, 0.35, 0.9)
e = (0.95, 0.15, 0.1)
ax2.quiver3D(e..., (s .- e)..., color = COLORS[5])

s = (0.55, 0.45, 0.1)
e = (0.75, 0.35, 0.9)
ax2.quiver3D(e..., (s .- e)..., color = COLORS[5])
ax2.text3D(1.1, 0.45, 0.1, "Stretch", color = COLORS[5])


fig.tight_layout()
fig.subplots_adjust(top = 0.95, wspace = 0.1)


# %% 01-test for chaos for standard map
using Statistics
sm = Systems.standardmap(;k = 1.0)
u1 = 0.1rand(2)
u2 = u1 .+ [π, 0]

fig, axs = subplots(2, 3;figsize = (figx, 2figy))
c = 0.2
for (i, u) in enumerate((u1, u2))
    x = trajectory(sm, 10000, u)[:, 2]
    axs[i, 1].plot(x[1:100]; color = "C$(i-1)", marker = "o", lw = 0.5)
    p, q = ChaosTools.trigonometric_decomposition(x, c)
    axs[i, 2].plot(p, q; color = "C$(i-1)", marker = "o", lw = 0.25, ms = 2)
    Z = ChaosTools.mmsd(mean(x), p, q, length(x)÷10, c)
    axs[i, 3].plot(Z; color = "C$(i-1)")
end
axs[1, 1].set_title("timeseries \$x_n\$")
axs[1, 2].set_title("process \$(p_n, q_n)\$")
axs[1, 3].set_title("sq. displ. \$Z_n\$")
fig.tight_layout(;pad = 0.25)

# %% 01-test chaos for henon heiles
using Statistics
sm = Systems.henonheiles()
u0s = (
    [0.0, -0.25, 0.42, 0.0], # chaotic
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0], # periodic
)

fig, axs = subplots(2, 3;figsize = (figx, 2figy))
c = 0.2
dt = 1.0
for (i, u) in enumerate(u0s)
    x = trajectory(sm, 10000, u; dt)[:, 2]
    axs[i, 1].plot(x[1:100]; color = "C$(i-1)", marker = "o", lw = 0.5)
    p, q = ChaosTools.trigonometric_decomposition(x, c)
    axs[i, 2].plot(p, q; color = "C$(i-1)", marker = "o", lw = 0.25, ms = 2)
    Z = ChaosTools.mmsd(mean(x), p, q, length(x)÷10, c)
    axs[i, 3].plot(Z; color = "C$(i-1)")
end
axs[1, 1].set_title("timeseries \$x_n\$")
axs[1, 2].set_title("process \$(p_n, q_n)\$")
axs[1, 3].set_title("sq. displ. \$Z_n\$")
fig.tight_layout(;pad = 0.25)
