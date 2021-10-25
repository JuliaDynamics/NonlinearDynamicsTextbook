# %% Recurrences
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, DynamicalBilliards, PyPlot
w = 0.2
pc = Particle(-0.01, 0.2, sqrt(3)/2)
pr = Particle(0.0, 1.9, 0.0)

ps = [pc, pr]
cs = ["C0", "C3"]
bd = billiard_mushroom(1, w, 1)

zbox = ((420, 130), (460, 170))

fig, axs = subplots(1, 3)
εs = [0.1, 0.1]
for i in 1:2
    ax = axs[i==1 ? 1 : 3]
    p = ps[i]
    N = 500
    bmap, = boundarymap(p, bd, N)
    tr = standardize(Dataset(bmap))
    R = RecurrenceMatrix(tr, εs[i])
    x, y = coordinates(R)
    ax.scatter(x, y, s = 1, color = cs[i])
    ax.set_aspect("equal")
    ℓ = dl_entropy(R)
    r = rt_entropy(R)
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)
    ax.set_xticks(0:250:500)
    ax.set_yticks(0:250:500)
    ax.grid()
    ax.set_title("\$H_\\ell = $(round(ℓ;digits= 1)),\\, H_r =$(round(r;digits=1))\$")
end
axs[3].set_yticklabels([])

# add zoomin
axis_zoomin!(axs[2], axs[1], zbox, zbox, "C1")
axs[2].axis("off")
axs[2].set_aspect("equal")
bmap, = boundarymap(pc, bd, 500)
tr = standardize(Dataset(bmap))
R = RecurrenceMatrix(tr, εs[1])
x, y = coordinates(R)
axs[2].scatter(x, y, s = 10, color = cs[1])
axs[2].set_xlim(zbox[1][1], zbox[2][1])
axs[2].set_ylim(zbox[1][2], zbox[2][2])

# Create indicative arrows
y0 = 137
yspan = 12
yc = (2y0 + yspan)/2
nice_arrow!(axs[2], 427, yc, 0, yspan; tex = "\$r=$(yspan)\$", color = "C2")


x0 = 445
y0 = 140
xspan = 10
yspan = 10
yc = (2y0 + yspan)/2
xc = (2x0 + xspan)/2
nice_arrow!(axs[2], xc, yc, xspan, yspan; tex = "\$\\ell=10\$", color = "C4")


add_identifiers!(fig, (axs[1], axs[3]))
fig.subplots_adjust(left = 0.08, bottom = 0.01, right = 0.95, top = 0.97, hspace = 0.1)
wsave(plotsdir("8","recurrence"), fig)
