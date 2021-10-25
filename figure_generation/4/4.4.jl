# Henon global bifurcation
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

he = Systems.henon(;a = 1.4, b = 0.3)

Z = 400
xg = range(-1.5, 1.5; length = Z)
yg = range(-0.5, 0.5; length = Z)

tr = trajectory(he, 10000; Ttr = 100)
ranges = (
    (range(-1.5, 1.5; length = Z), range(-0.5, 0.5; length = Z)),
    (range(-1.25, -0.85; length = Z), range(-0.45, -0.3; length = Z)),
)

# calculate basins
allbasins = []
for (i, (xg, yg)) ∈ enumerate(ranges)
    i > 1 && break
    @time basins, attractors = basins_of_attraction((xg, yg), he; 
    mx_chk_fnd_att = 30, mx_chk_lost = 500, horizon_limit = 1e2)
    push!(allbasins, basins)
end
push!(allbasins, allbasins[1])

# %% Plot
fig, axs = subplots(1, 3)

for (i, (xg, yg)) ∈ enumerate(ranges)
    basins = allbasins[i]
    LC =  matplotlib.colors.ListedColormap
    cmap = LC([matplotlib.colors.to_rgb(c) for c in ("C1", (1.0, 1.0, 0.9))])
    vmin = -1; vmax = 1
    # axs[i].pcolormesh(xg, yg, basins'; cmap, vmin, vmax)
    axs[i].pcolormesh(ranges[1][1], ranges[1][2], basins'; cmap, vmin, vmax)
    # get a quality attractor to plot
    axs[i].plot(columns(tr)...; marker = "o", ms = 4, color = "C2", ls = "None")
    axs[i].set_xlim(xg[1], xg[end])
    axs[i].set_ylim(yg[1], yg[end])
    if i == 1
        axs[i].set_xlabel("\$x\$", labelpad=-20)
        axs[i].set_ylabel("\$y\$", labelpad=-40)
        axs[i].set_xticks([xg[1], xg[end]])
        axs[i].set_yticks([yg[1], yg[end]])
    end
end
axs[2].axis("off")
# add rectangles
xg, yg = ranges[2]
for i in 1:2
    line, = axs[i].plot(
        [xg[1], xg[end], xg[end], xg[1], xg[1]],
        [yg[1], yg[1], yg[end], yg[end], yg[1]],
        color="C3", lw = 3
    )
    line.set_clip_on(false)
end

# Transient chaos trajectory
he = Systems.henon([0.1, 0.1]; a = 1.427, b = 0.3)
N = 2700
tr = trajectory(he, N, Ttr = 10)
for i in 1:2
    axs[i].plot(columns(tr)...; marker = "o", ms = 2, ls = "None", color = "C0")
end

axs[3].plot(0:N, tr[:, 1], marker = "o", ls = "None", ms = 5, color = "C0")
axs[3].set_ylim(-4, 2)
axs[3].set_xticks([0, 2500])
axs[3].set_xlabel("\$n\$", labelpad = -20)
axs[3].set_ylabel("\$x_n\$"; labelpad = -20)
add_identifiers!(fig)
fig.tight_layout(pad=0.3)

wsave(plotsdir("4", "boundary_crisis"), fig)