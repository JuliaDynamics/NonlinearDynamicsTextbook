# Period doubling bifurcation + cobweb for logistic map
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

lo = Systems.logistic(0.4)
rs = (2.8, 3.3)

fig, axs = subplots(1, 2, figsize = (0.66figx, figy))

T = 1000
xs = 0:0.0001:1

function cobweb(t) # transform timeseries t into cobweb (point2D)
    cx = Float64[]; cy = Float64[]
    for i ∈ 1:length(t)-1
        push!(cx, t[i]); push!(cy, t[i])
		push!(cx, t[i]); push!(cy, t[i+1])
    end
    # add line to x axis
    pushfirst!(cx, cx[1])
    pushfirst!(cy, 0)
    return cx, cy
end

for (i, r) in enumerate(rs)
    ax = axs[i]
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xticks([0, 1])
    ax.set_xlabel("\$x\$"; labelpad=-20)
	set_parameter!(lo, 1, r)
	f  = lo.f.(xs, Ref([r]), 0)
    f² = lo.f.(f, Ref([r]), 0)
	ax.plot(xs, f, color = "C0")
	ax.plot([0, 1], [0,1]; lw = 2, color = "C3", ls = "--")
	x = trajectory(lo, T)
    cx, cy = cobweb(x)
    if i>1
        ax.plot(xs, f², color = "C1")
        ax.text(0.15, 0.875, "\$f^{(2)}\$"; color = "C1", transform = ax.transAxes, size = 32)
    end
    ax.plot(cx, cy; lw = 1.0, color = "C2")
    ax.text(0.05, 0.875, "\$f\$"; color = "C0", transform = ax.transAxes, size = 32)
end
axs[2].set_yticklabels([])

fig.tight_layout(pad = 0.37)
add_identifiers!(fig)

# wsave(plotsdir("period_doubling_bifurcation"), fig)