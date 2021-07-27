# %% CCM demonstration
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random
using Neighborhood, Statistics, LinearAlgebra

function ccm(x, y, d, τ, w = τ)
    Mx = embed(x, d, τ); theiler = Theiler(w); tree = KDTree(Mx)
    idxs = bulkisearch(tree, Mx, NeighborNumber(d+1), theiler)
    ỹ = copy(y)
    for i in 1:length(Mx)
        J = idxs[i]
        xᵢ = Mx[i]; n1 = norm(xᵢ - Mx[J[1]])
        w = [exp(-norm(xᵢ - Mx[j])/n1) for j in J]
        w ./= sum(w)
        ỹ[i] = sum(w[k]*y[j] for (k, j) in enumerate(J))
    end
    return cor(y, ỹ)
end

ds = Systems.lorenz()
tr = trajectory(ds, 1000; dt = 0.1, Ttr = 100)
x, y, z = columns(tr)
τx = estimate_delay(x, "mi_min")
τy = estimate_delay(y, "mi_min")

Ns = 50:50:3000
ρx = [ccm(x[1:N], y[1:N], 3, 5) for N in Ns]
ρy = [ccm(y[1:N], x[1:N], 3, 5) for N in Ns]
ρz = [ccm(z[1:N], y[1:N], 3, 5) for N in Ns]

# %%

fig = figure(figsize = (figx/2, figy))
ax = gca()
ax.plot(Ns, ρx; label = "\$x \\to y\$")
ax.plot(Ns, ρy; label = "\$y \\to x\$", color = "C2")
ax.set_ylim(0.68, 1)
ax.set_xticks(0:1000:3000)
ax.legend(loc = "lower right")
ax.set_ylabel("CCM measure")
ax.set_xlabel("\$N_s\$"; labelpad = -15)

fig.tight_layout(pad=0.3)
wsave(plotsdir("ccm_result"), fig)
