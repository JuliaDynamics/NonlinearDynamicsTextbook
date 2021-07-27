# Henon fractal zoom
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

fig, axes = subplots(1,3, figsize = (figx,figx/4))
botrow = axes

he = Systems.henon()
N = 10000
tr = trajectory(he, N, Ttr = 100)
integ = integrator(he)
data1 = Dataset([integ.u])
data2 = Dataset([integ.u])
n, m = 1, 1

zbox1 = ((0.93, 0.05), (1.17, 0.13))
zbox2 = ((1.1275, 0.078), (1.1425, 0.085))
while n < N
    step!(integ)
    if zbox1[1][1] < integ.u[1] < zbox1[2][1] && zbox1[1][2] < integ.u[2] < zbox1[2][2] && m < N
        push!(data1, integ.u)
        global m+=1
    end
    if zbox2[1][1] < integ.u[1] < zbox2[2][1] && zbox2[1][2] < integ.u[2] < zbox2[2][2]
        push!(data2, integ.u)
        global n+=1
    end
end

for ax in botrow[1:2]
    ax.clear()
    ax.plot(tr[:, 1], tr[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1)
end
botrow[1].set_xticks([])
botrow[1].set_yticks([])
botrow[2].clear()
botrow[2].plot(data1[:, 1], data1[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1)
botrow[2].axis("off")
botrow[3].clear()
botrow[3].plot(data2[:, 1], data2[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1)
botrow[3].axis("off")

for (i, zbox) in enumerate((zbox1, zbox2))
    axis_zoomin!(botrow[i+1], botrow[i], zbox, zbox, COLORS[i+1])
    botrow[i+1].set_xlim(zbox[1][1], zbox[2][1])
    botrow[i+1].set_ylim(zbox[1][2], zbox[2][2])
end

fig.tight_layout(pad = 0.3)
fig.subplots_adjust(top = 0.98, bottom = 0.05, right = 0.95, left = 0.05, wspace=0.1, hspace = 0.1)
wsave(plotsdir("henon_zoom"), fig)
