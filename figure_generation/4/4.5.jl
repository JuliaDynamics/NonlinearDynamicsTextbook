# period doubling route
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

fig = figure(figsize = (figx, 2figy))
axlo = fig.add_subplot(3, 1, 1)
axro = fig.add_subplot(3, 1, 3)

zoomaxs = []
for i in 1:3
	push!(zoomaxs, fig.add_subplot(3, 3, 3+i))
end

axlo.clear()
ds = Systems.logistic(0.5)
i = 1
pvalues = range(2.5, 4; length = 2001)
n = 2000
Ttr = 2000
p_index = 1
output = orbitdiagram(ds, i, p_index, pvalues; n = n, Ttr = Ttr)

L = length(pvalues)
x = Vector{Float64}(undef, n*L)
y = copy(x)
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= pvalues[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end
axlo.plot(x, y, ls = "None", ms = 0.5, color = "black", marker = "o", alpha = 0.01)
axlo.set_xlim(pvalues[1], pvalues[end]); axlo.set_ylim(0,1)
axlo.set_ylabel("\$x\$", labelpad = 10)
axlo.set_xticks(2.5:0.5:4)
axlo.set_yticks(0:0.2:1)
axlo.set_xlabel("\$r\$"; labelpad = -20)
axlo.text(2.6, 0.1, "logistic map", size = 36)
r∞ = 3.56995
axlo.axvline(r∞; color = "C4", alpha = 0.7, ls = "dashed")

# Zoomin boxes
zbox1 = ((3.629, 0.295), (3.640, 0.32))
zbox2 = ((3.6338, 0.304), (3.6342, 0.3065))
zbox3 = ((3.63407, 0.30495), (3.63408, 0.3052))
allboxes = (zbox1, zbox2, zbox3)
for (k, zbox) in enumerate(allboxes)
	pval = range(zbox[1][1], zbox[2][1]; length = 1001)
	nsm = n÷3
	output = orbitdiagram(ds, i, p_index, pval;
		n=nsm, Ttr, ulims = (zbox[1][2], zbox[2][2])
	)
	L = length(pval)
	x = Vector{Float64}(undef, nsm*L)
	y = copy(x)
	for j in 1:L
	    x[(1 + (j-1)*nsm):j*nsm] .= pval[j]
	    y[(1 + (j-1)*nsm):j*nsm] .= output[j]
	end
	zoomaxs[k].clear()
	zoomaxs[k].plot(x, y, ls = "None", ms = 0.5, color = "black", marker = "o", alpha = 0.01)
	zoomaxs[k].set_xlim(zbox[1][1], zbox[2][1])
	zoomaxs[k].set_ylim(zbox[1][2], zbox[2][2])
	zoomaxs[k].axis("off")
end

zbox0 = ((pvalues[1], 0.0), (pvalues[end], 1.0))
allzoomaxs = [axlo, zoomaxs...]
for i in 1:3
	origin = allzoomaxs[i]
	zoomin = allzoomaxs[i+1]
    zbox = allboxes[i]
	c = ("C1", "C3", "C5")[i]
    axis_zoomin!(zoomin, origin, zbox, zbox, c;	connect_lines = i > 1, lw = 4.0)
end

# arrow indicator
axlo.arrow(
	(zbox1[1][1] + zbox1[2][1])/2, 0.05, 0, zbox1[2][2] - 0.1;
	color = "C1", width = 0.005, length_includes_head = true, head_length = 0.1
)

# Roessler
ro = Systems.roessler()
pvalues = range(1.5, stop = 5.18, length = 2001)
i = 1
plane = (2, 0.0)
tf = 4000.0
p_index = 3
output = produce_orbitdiagram(ro, plane, i, p_index, pvalues;
                              tfinal = tf, Ttr = 2000.0)

for (j, p) in enumerate(pvalues)
    axro.plot(fill(p, length(output[j])), output[j], lw = 0,
    marker = "o", ms = 0.2, color = "black", alpha = 0.1)
end
axro.set_xlim(pvalues[1], pvalues[end]);
axro.set_ylim(minimum(minimum(o) for o in output), maximum(maximum(o) for o in output))
axro.set_ylabel("\$x\$ (at section)", labelpad = 25)
axro.set_xlabel("\$c\$"; labelpad = -20)
axro.set_yticks(3:3:9)
axro.text(1.7, 7, "Rössler system", size = 36)
fig.tight_layout(pad=0.3)
fig.subplots_adjust(hspace = 0.3, wspace = 0.1)
add_identifiers!(fig, (axlo, axro))
wsave(plotsdir("4", "orbit_diagrams"), fig)