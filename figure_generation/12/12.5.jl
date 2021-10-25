# %% Carbon cycle extinction events
# model of Rothman, 10.1073/pnas.1905164116
# main dynamic rule, eqs.(11, 12) of paper
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Statistics, DelimitedFiles
function excitable_carbon(u, p, t)
    c, w = u
    μ, b, θ, cₓ, cₚ, ν, w₀, γ, f₀, c_f, β = p
    s = excitable_carbon_s(c, cₚ, γ)
    sbar = 1 - excitable_carbon_s(c, cₓ, γ)
    f = excitable_carbon_f(c, f₀, β, c_f)
    cdot = μ*(1 - b*s - θ*sbar - ν) + w - w₀
    cdot *= f
    wdot = μ*(1 - b*s + θ*sbar + ν) - w + w₀
    return SVector(cdot, wdot)
end
# functions used in equations of motion
excitable_carbon_f(c, f₀, β, c_f) = f₀*(c^β / (c^β + c_f^β))
excitable_carbon_s(c, cₚ, γ) = c^γ/(c^γ + cₚ^γ)

u0 = [100.0, 2800.0] # random, from Fig. 4

# Parameter default values in Table 1 of appendix, but cₓ same as Fig. 5
p0 = [250, 4, 5, 55, 110, 0, 2000, 4, 0.694, 43.9, 1.7]
# μ, b, θ, cₓ, cₚ, ν, w₀, γ, f₀, c_f, β = p0


ec = ContinuousDynamicalSystem(excitable_carbon, u0, p0)

# streamplot:
xgrid = 0:5:190
ygrid = 2000:20:3500
ux = zeros(length(xgrid), length(ygrid))
uy = copy(ux)
for (i, x) in enumerate(xgrid)
    for (j, y) in enumerate(ygrid)
        ux[i, j], uy[i, j] = ec.f(SVector(x, y), ec.p, 0)
    end
end

fig, axs = subplots(1,2; figsize = (0.66figx, figy))

axs[1].streamplot(Vector(xgrid), Vector(ygrid) ./1000 , ux', uy' ./ 1000;
    linewidth = 1.0, density = 0.25, color = "C4", arrowsize = 1.0
)

# plot trajectories
c_stable = 83.58 # approximate value
for (i, u) in enumerate(([c_stable, 2075.0], [c_stable, 2102.0], ))
    tr = trajectory(ec, 100.0, u)
    c, w = columns(tr)
    w ./= 1000
    for ax in axs
        ax.plot(c, w; color = "C$(i-1)")
    end
end
ax = axs[1]
ax.set_xlabel("\$c\$"; labelpad = -15)
ax.set_ylabel("\$w \\times 10^{-3}\$")
ax.set_xticks(0:60:180)

zbox = ((65, 2.05), (100, 2.25))

# Plot fixed point and arrows
fp = trajectory(ec, 10.0, [83.56, 2075]; Ttr = 1000)[end]
axs[2].plot([fp[1]], [fp[2]/1000]; marker = "o", mew=4, ms = 15, mec="C2", zorder = 99, mfc="white")


# axs[2].arrow(
#     fp[1], fp[2]/1000, 0,  2.1003 - fp[2]/1000,
#     zorder = 99, width = 0.01,
# )

# axs[2].arrow(
#     fp[1], 2.102, 0, fp[2]/1000 - 2.102;
#     zorder = 99, width = 0.01, length_includes_head = true,
#     # head_width = 2,
#     head_starts_at_zero = true
# )



axis_zoomin!(axs[2], axs[1], zbox, zbox, "C3"; α = 0.75)
axs[2].axis("off")
fig.tight_layout(pad = 0.3)
fig.subplots_adjust(wspace = 0.1)

# add_identifiers!(fig)
wsave(plotsdir("12", "carboncycle"), fig)
