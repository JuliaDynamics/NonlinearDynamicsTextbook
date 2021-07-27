# Intermittency simple Manneville map
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot
lo = Systems.manneville_simple(0.454)

rs = [-0.1, 0.0001, 0.05]

fig = figure(figsize = (figx, 1.5figy))
ax1 = fig.add_subplot(3, 2, 1)
ax2 = fig.add_subplot(3, 2, 3; sharex = ax1)
ax3 = fig.add_subplot(3, 2, 5; sharex = ax1)
axs = (ax1, ax2, ax3)
axc = fig.add_subplot(1, 2, 2)
add_identifiers!(fig)

T = 180
xs = 0:0.0001:1

for i in 1:length(rs)
	set_parameter!(lo, 1, rs[i])
	x = trajectory(lo, T; Ttr = 100)
	chaotici = findall(a -> a>0.1, x)
	regulari = setdiff(1:length(x), chaotici)
	axs[i].plot(0:length(x)-1, x; lw = 0.5, color = "C0")
	axs[i].plot(regulari .- 1, x[regulari]; ls = "None", marker = "o", color = "C0")
	axs[i].plot(chaotici .- 1, x[chaotici]; ls = "None", marker = "o", color = "C1")
	axs[i].set_ylim(-0.1, 1.2)
	axs[i].set_yticks([0, 1])
	# axs[i].text(0.99, 0.90, "\$r=$(rs[i])\$"; bbox = bbox,
	# transform = axs[i].transAxes, va = :top, ha=:right, fontsize = 22)
end
# axs[1].set_xlim(0, T)
for ax in (ax1, ax2); ax.tick_params(labelbottom=false); end
axs[1].set_xlim(0, T)
axs[1].set_xticks(0:60:180)
axs[3].set_xlabel("\$n\$"; labelpad = -25)
axs[2].set_ylabel("\$x_n\$"; labelpad = -10)
axs[2].axvspan(21, 40; color = "C3", alpha = 0.5)

# laminar arrow
xc = (135 + 169)/2
xspan = (169 - 135)
nice_arrow!(ax2, xc, 1.0, xspan, 0; tex = "\$\\ell\$", yo = 0.0, xo = -xspan/2 - 5)

xc = (65 + 100)/2
xspan = abs(65 - 100)
nice_arrow!(ax2, xc, 1.0, xspan, 0; tex = "\$c\$", yo = 0.0, xo = xspan/2, color = "C1")

# Do the cobweb plot
function cobweb(t) # transform timeseries t into cobweb (point2D)
    cx = Float64[]; cy = Float64[]
    for i ∈ 1:length(t)-1
        push!(cx, t[i]); push!(cy, t[i])
		push!(cx, t[i]); push!(cy, t[i+1])
    end
    return cx, cy
end


axc.set_xlim(0.0, 1.0)
axc.set_ylim(0.0, 1.0)
axc.set_xticks([0,1])
axc.set_yticks([0,1])
axc.set_ylabel("\$x_{n+1}\$"; labelpad = -15)
axc.set_xlabel("\$x_{n}\$"; labelpad = -25)
axi = axc.inset_axes([0.05, 0.55, 0.3, 0.4])

for ax in (axc, axi)
	_r = rs[2]
	set_parameter!(lo, 1, _r)
	f = lo.f.(xs, Ref([_r]), 0)
	ax.plot(xs, f, color = "C2")
	ax.plot([0, 1], [0,1]; lw = 1.0, color = "k")
	x = trajectory(lo, 60; Ttr = 100)
	cx, cy = cobweb(x)
	ax.plot(cx, cy; lw = 1.0, alpha = 0.5, color = COLORS[1])
end

z = 0.01
w = 0.05
zbox = ((z, z), (z+w, z+w))
axis_zoomin!(axi, axc, zbox, zbox, "C3"; dir = :top)
axi.set_ylim(z, z+w)
axi.set_xlim(z, z+w)
axi.axis("off")

fig.tight_layout(pad = 0.3)
# wsave(plotsdir("intermittency_manneville"), fig)
