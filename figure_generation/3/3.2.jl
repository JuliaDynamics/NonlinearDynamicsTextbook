# definition of lyapunov
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

using LinearAlgebra
lo = Systems.lorenz([20,20,20.0])
X₁ = trajectory(lo, 45)
u₂ = get_state(lo) .+ 1e-6
X₂ = trajectory(lo, 45, u₂)
δ  = norm.(X₂.data .- X₁.data)
λ = lyapunov(lo, 100000)

fig = figure(figsize=(figx/2,figy/2))
ax = gca()
ax.plot(0:0.01:45, log.(δ), c = COLORS[1], label ="\$\\ln(\\delta(t)))\$")
ax.set_yticks(-12:4:4)
ax.set_xlabel("\$t\$", labelpad=-20)
# Lyapunov
ax.plot([0, 15] .+ 4, λ .* [0, 15] .- 13, color = COLORS[2])
ax.text(12, -9, "\$\\lambda_1\$=$(round(λ;digits=2))", color = COLORS[2])
ax.legend(;handlelength = 1.0)
xticks(0:15:45)
tight_layout(pad = 0.25)
subplots_adjust(bottom = 0.24, left = 0.12)
wsave(plotsdir("lyapunov"), fig)
