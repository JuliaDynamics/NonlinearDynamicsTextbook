# Koch snowflake
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

linepoints = SVector{2}.([[0.0; 0.0], [1.0; 0.0]])
flakepoints = SVector{2}.([[0.0; 0.0], [0.5; sqrt(3)/2], [1; 0.0], [0.0; 0.0]])
function koch(points, maxk, α = sqrt(3)/2)
  Q = @SMatrix [0 -1; 1 0]
  for k = 1:maxk
    n = length(points)
    new_points = eltype(points)[]
    for i = 1:n-1
      p1, p2 = points[i], points[i+1]
      v = (p2 - p1) / 3
      q1 = p1 + v
      q2 = p1 + 1.5v + α * Q * v
      q3 = q1 + v
      push!(new_points, p1, q1, q2, q3)
    end
    push!(new_points, points[end])
    points = new_points
  end
  return points
end

# Plot construction
fig, axes = subplots(2,3, figsize = (figx, 2figy))
toprow = axes[1:2:5]
botrow = axes[2:2:6]
kochcolors = COLORS[1:4]

for (o, ax) in enumerate(toprow)
    ax.clear()
    ax.set_aspect("equal")
    k1 = koch(flakepoints, o-1)
    k2 = koch(flakepoints, o)

    ax.plot([a[1] for a in k1], 0.5 .+ [a[2] for a in k1], kochcolors[o], lw = 2.5)
    ax.plot([a[1] for a in k2], 0.5 .+ [a[2] for a in k2], kochcolors[o+1], 
    label="Step $(o)", alpha = 0.75, lw = 3)
    # ax.legend(loc = "upper left", framealpha = 1.0, ncol = 2, handlelength = 1)
    ax.text(-0.1, 0.88, "step $(o)", transform = ax.transAxes, size=28)
    # ax.set_ylim(ys...)
    # ax.set_xlim(xs...)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis("off")
end

# Plot versus circle
points1 = [[0.5; 0.865], [2/3, 0.577361]]
points1 = [[0.5; 0.577361], [2/3, 0.865]]
points2 = [[0.537233, 0.801228], [0.555887, 0.768912]]
points2 = [[0.537233, 0.768912], [0.555887,0.801228]]
allpoints = [flakepoints, points1, points2]

for i in 1:3
    botrow[i].clear()
    largekoch = koch(allpoints[i], 6)
    kx, ky = [a[1] for a in largekoch], [a[2] for a in largekoch]
    botrow[i].plot(kx, ky, lw = 1.0)
    botrow[i].set_aspect("equal")
    botrow[i].add_artist(plt.Circle(
        (0.5,0.28900), 0.5, lw=2.0, color="C1", alpha = 0.75, fill=false
    ))
    botrow[i].axis("off")
end
for i in 1:2
    origin = botrow[i]
    zoomin = botrow[i+1]
    zbox = allpoints[i+1]
    axis_zoomin!(zoomin, origin, zbox, zbox, "C$(i+1)"; α = 0.5)
end

fig.tight_layout(pad = 0.1)
fig.subplots_adjust(wspace = 0.01, hspace = 0.01, bottom = 0.01, left = 0.01)
# wsave(plotsdir("koch"), fig)
