# Magnetic pendulum basins of attraction
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

α = 0.2; ω = 1.0; d = 0.3
gx = gy = range(-5, 5; length = 1500)
config = @strdict(α, ω, d, gx, gy)

function magnetic_basins(config)
    @unpack α, ω, d, gx, gy = config
    ma = Systems.magnetic_pendulum(; α, ω, d)
    if haskey(config, "attractors")
        attractors = config["attractors"]
        @time basins, attractors = basins_of_attraction((gx, gy), ma;
        mx_chk_att = 1, mx_chk_fnd_att = 10, attractors)
    else
        @time basins, attractors = basins_of_attraction((gx, gy), ma;
        mx_chk_att = 1, mx_chk_fnd_att = 10, mx_chk_lost = 10)
    end
    return @strdict(gx, gy, basins, attractors)
end

# produce high quality
file, _ = produce_or_load(datadir(), config, magnetic_basins; prefix = "magnetic_basins", tag = false)
attractors = file["attractors"]

# produce zoomed version
gx = range(1.80, 1.95; length = 1000)
gy = range(0, 0.12; length = 1000)
config = @strdict(α, ω, d, gx, gy, attractors)
produce_or_load(datadir(), config, magnetic_basins; prefix = "magnetic_basins_zoomed", tag = false)

# %% Plot this
@unpack gx, gy, basins, attractors = wload(datadir(savename("magnetic_basins", config, "jld2")))
ma = Systems.magnetic_pendulum(; α, ω, d)

LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(0.33figx, figy))
ax = gca()
ax.pcolormesh(gx, gy, basins'; cmap, shading = "gouraud")
ax.set_aspect("auto")
xticks([-5, 5])
yticks([-5, 5])
xlabel("\$x\$", labelpad=-32)
ylabel("\$y\$", labelpad=-32)
for m in values(attractors)
    m1, m2 = columns(m)
    scatter(m1, m2; color = "white", edgecolor = "black", zorder = 99, s = 100)
end
# wsave(plotsdir("magneticpendulum"), fig)

# %% Plot zoomed version
@unpack gx, gy, basins, attractors = wload(datadir(savename("magnetic_basins_zoomed", config, "jld2")))
# fig = figure(figsize=(figx/2, figx/2))
x0 = 0.55; wx = 0.42; y0 = 0.65; wy = 0.32;
axin = ax.inset_axes([x0, y0, wx, wy])

axin.pcolormesh(gx, gy, basins'; cmap, shading = "gouraud")
# axin.set_aspect("equal")

# xticks([1.8, 1.95], size = 20)
# yticks([0, 0.12], size = 20)
zbox = ((gx[1], gy[1]), (gx[end], gy[end]))
axis_zoomin!(axin, ax, zbox, zbox, "C3"; dir = :top, lw = 4)

axin.axis("off")

fig.tight_layout(pad=0.3)

wsave(plotsdir("5", "magneticpendulum.png"), fig)