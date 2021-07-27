# %% Demonstration of delay embeddings
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

lo = Systems.lorenz([0, 10, 0.0])
tr = trajectory(lo, 1000; Ttr=10, reltol = 1e-12, abstol = 1e-12)
x, y, z = columns(tr)
w = x
τ = estimate_delay(w, "mi_min")
R = embed(w, 3, τ)

close("all")
N = 5000
fig = figure()
ax1 = subplot(131, projection = "3d")
ax1.plot3D(x[1:N], y[1:N], z[1:N], lw = 1.0, color = "C0")
ax1.set_title("original set \$A\$", pad = 30)
ax = ax1
for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end
# ax1.text(-20, -20, 0, "\$D_0(A)\$ = $(round(D_A;digits=3))", color = "C1")

ax2 = subplot(132)
wx = 2400:3200
ax2.plot(wx, w[wx])
midpoint = (wx[1] + wx[end])/2
tau = midpoint:(midpoint+τ)
ax2.plot(tau, fill(10, length(tau)), color = "C3")
ax2.text(midpoint, 7.5, "\$\\tau\$", color = "C3", size = 32)
ax2.set_title("measurement \$w=x\$")
# make fancy axis
ax2.spines["bottom"].set_position("center")
ax2.spines["right"].set_color("none")
ax2.spines["top"].set_color("none")
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel("\$t\$")
ax2.set_ylabel("\$w\$", rotation = 0)
ax2.xaxis.set_label_coords(1.0, 0.5)
ax2.yaxis.set_label_coords(-0.05, 0.9)

ax3 = subplot(133, projection = "3d")
ax3.plot(R[1:N, 1], R[1:N, 2], R[1:N, 3], color = "C1", lw = 1.0)
ax = ax3
for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end
ax3.set_title("reconstruction \$R\$", pad = 30)
# ax3.set_title("reconstruction \$R\$\n\$(\\gamma=$(γ),\\tau=$(τ))\$")
# ax3.text(-5,-25,-10, "\$D_0(R)\$ = $(round(D_R;digits=3))", color = "C2")
ax1.dist = 9
ax3.dist = 9

fig.tight_layout(pad = 0.3)
# wsave(plotsdir("delayembedding"), fig)