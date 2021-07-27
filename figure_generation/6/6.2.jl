# %% value of τ demonstration
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

d = 3
τ1 = 1
τ2 = 17
τ3 = 3τ2
θ = -30
e = 10
lo = Systems.lorenz([0, 10, 0.0])
A = trajectory(lo, 1000; Ttr=1000, reltol = 1e-12, abstol = 1e-12)
w = A[:, 1]
close("all")

fig = figure()
for (i, τ) in enumerate((0, τ1, τ2, τ3))
    R = τ == 0 ? A : embed(w, d, τ)
    ax = subplot(1, 4, i, projection = "3d")
    ax.plot(R[1:N, 1], R[1:N, 2], R[1:N, 3],
    color = i == 1 ? "C0" : "C1", lw = 1,
    zorder = 1)
    ax.view_init(elev=e, azim=θ)
    ax.set_title(i == 1 ? "original" : "\$R:\\,\\tau=$(τ)\$", pad = -10)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end
    # dj = 167
    # for j in (7827, )
    #     χ, ψ, ζ = R[j:dj:j+dj, 1], R[j:dj:j+dj, 2], R[j:dj:j+dj, 3]
    #     ax.scatter3D(χ, ψ, ζ,depthshade = false,c = "C2", zorder = 99, s = 50)
    # end
    ax.dist = 9
end

# ax = subplot(1, 4, 1, projection = "3d")
# ax.plot3D(x[1:N], y[1:N], z[1:N], lw = 1.0, color = "C1")
# ax.view_init(elev=e, azim=θ)
# for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end

fig.tight_layout()
fig.subplots_adjust(wspace = 0.0001, bottom = -0.1, top = 1, left = 0.0, right = 1)
# wsave(plotsdir("taudemonstration"), fig)
