# %% Chaotic scattering
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using PyPlot, DynamicalBilliards
fig, axs = subplots(1,3)
axs[3].axis("off")
ax3 = fig.add_subplot(1, 6, 5)
ax4 = fig.add_subplot(1, 6, 6)
DynamicalBilliards.obcolor(::Obstacle) = matplotlib.colors.to_rgb(COLORS[3])

# First plot: illustration of scattering function
# I GIVE UP, I'll Make this plot in PowerPoint
off = 0.3
r = 0.2
offset = [0.0, off]
center = 2.5
R(φ) = [cos(φ) -sin(φ);
        sin(φ)  cos(φ)]
R3 = R(2π/3)
enclosing = billiard_sinai(r, 2center, 2center)
enclosing = enclosing[2:5]

# disk0 = Disk([center,0], 1.5r, "0")
# sca(axs[1])
# plot(disk0)
# axs[1].set_xlim(-1.5 + center, 0.5 + center)
# axs[1].set_ylim(-1+0.3, 1+0.3)
axs[1].axis("off")

# p0 = Particle(center-1, 0.22, 0)
# bd0 = Billiard(enclosing..., disk0)
# x,y,vx,vy = timeseries(p0, bd0, 2; dt = 0.01)
# axs[1].plot(x[1:130],y[1:130]; color = "C1")
# x0, y0 = p0.pos
# vx0, vy0 = p0.vel
# axs[1].quiver(x0, y0, 0.08vx0, 0.08vy0; angles = "xy", scale = 0.5, width = 0.018, color="C1", zorder = 99)
# axs[1].text(x0, y0-0.2, "\$\\phi\$")
# axs[1].text(x0, y0-0.2, "\$a\$")

# Okay now clarify the scattering function
disk1 = Disk([center, center] .+ offset, r, "green")
disk2 = Disk([center, center] .+ R3*offset, r, "red")
disk3 = Disk([center, center] .+ R3*R3*offset, r, "purple")

disks = Billiard(disk1, disk2, disk3)
plot(disks; ax = axs[2])
axs[2].set_xlim(1.9,3.1)
axs[2].set_ylim(2.0,3.2)
scattering = Billiard(enclosing..., disk1, disk2, disk3)

axs[2].axis("off")

# make some particles

terminate(n, τ, i, p) = i ∈ 1:4
φs = 5π/6 .+ range(0; step = 0.001, length = 3) .+ 0.17

ps = [Particle((R(φ)*[0.5, 0] .+ [center, center])..., φ+π) for φ in φs]

partcol = ( "C3", "C1", "C0",)
for (i, p) in enumerate(ps);
    x0, y0 = p.pos
    vx0, vy0 = p.vel
    axs[2].quiver(x0, y0, 0.04vx0, 0.04vy0; angles = "xy", scale = 0.5, width = 0.01, 
    color = partcol[i], zorder = 99)
    x,y = timeseries(p, scattering, terminate)
    axs[2].plot(x, y; color = partcol[i], lw = 1, alpha = 0.75)
end


# Detailed computation of input-output function
φs = range(0, 2π/3; length = 100_000)
θs = zero(φs)

for (i, φ) in enumerate(φs);
    p = Particle((R(φ)*[0.5, 0] .+ [center, center])..., φ+π)
    timeseries!(p, scattering, terminate)
    θs[i] = atan(p.vel[2], p.vel[1])
end

ax3.plot(φs, θs .+ π; ls = "None", marker = "o", ms = 0.5, alpha = 0.5)
ax3.set_xticks([0, 2π/3])
ax3.set_xticklabels(["0", "2π/3"])
ax3.set_yticks([0, 2π])
ax3.set_yticklabels(["0", "2π"])
ax4.set_yticks([])
ax3.set_xlim(0, 2π/3)
ax4.set_ylim(0, 2π)
ax3.set_ylim(0, 2π)

φs = range(0.1555, 0.1567; length = 100_000)
θs = zero(φs)
for (i, φ) in enumerate(φs);
    p = Particle((R(φ)*[0.5, 0] .+ [center, center])..., φ+π)
    timeseries!(p, scattering, terminate)
    θs[i] = atan(p.vel[2], p.vel[1])
end
ax4.plot(φs, θs .+ π; ls = "None", marker = "o", ms = 0.5, alpha = 0.5)
ax4.set_xlim(φs[1], φs[end])
ax4.set_xticks([])
ax3.set_xlabel("\$\\phi_i\$", labelpad = -25)
ax3.set_ylabel("\$\\phi_o = \\mathbf{S}(\\phi_i)\$", labelpad = -25)
# ax3.axvspan(φs[1], φs[end]; color = "C1") # not even visible lol

fig.tight_layout(;pad = 0.25)
fig.subplots_adjust(wspace = 0.1)

# wsave(plotsdir("chaoticscattering"), fig)