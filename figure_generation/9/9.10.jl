using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random, OrdinaryDiffEq

a = 0.2
b = 0.2
c = 5.7
amu = 0.02
ω1 = 1. + amu
ω2 = 1. - amu
k1 = 0.01
k2 = k1

y_init = [0.0, 0.2, 0., 0.11, 0.19, 0.1]

croesslers = Systems.coupled_roessler(y_init; ω1, ω2, a, b, c, k1, k2)
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9, maxiters = Inf)
Δt = 0.01
T = 1000
t = 0:Δt:T
L = length(t)÷10

ks = [0.03, 0.04, 0.05] # estimated from the Lyapunov exponents plot

bbox2 = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 0.9)
axs_for_labels = []

fig = figure(figsize = (figx, 1.5figy))
for (i, k) in enumerate(ks)
    M = length(ks)
    @show k 
    set_parameter!(croesslers, 6, k)
    set_parameter!(croesslers, 7, k)
    X = trajectory(croesslers, T; Ttr = 5000, Δt, diffeq...)

    row = (i-1)*3
    ax1 = fig.add_subplot(M, 3, (1, 2) .+ row)
    ax2 = fig.add_subplot(M, 3, 3+row)

    ax1.plot(t[1:L], X[1:L, 1])
    ax1.plot(t[1:L], X[1:L, 4]; lw = 3.0, color = COLORS[2])
    ax2.plot(X[:, 1], X[:, 4]; lw = 1.0, color = COLORS[1])
    push!(axs_for_labels, ax1)

    λs = lyapunovspectrum(croesslers, 50000, 4; Ttr = 1000, diffeq...)
    # ax1.text(0.02, 0.9, "\$k = $(k)\$\n\$ \\lambda_1 = $(round(λs[1];digits=2))\$\n\$\\lambda_2 = $(round(λs[2];digits=2))\$",
    #     transform = ax1.transAxes, va = :top, bbox = bbox2, fontsize = 28
    # )
    ax1.text(0.02, 0.5, "\$k = $(k),\\, \\lambda_1 = $(round(λs[1];digits=2))\\,\\lambda_2 = $(round(λs[2];digits=2))\$",
        transform = ax1.transAxes, va = :center, bbox = bbox2, fontsize = 28
    )
    @show λs[1:4]
    ax1.set_xlim(0, t[1:L][end])
    ax1.set_ylim(-12,12)
    ax2.set_xlim(-12,12)
    ax2.set_ylim(-12,12)
    ax2.set_yticklabels([])
    ax1.set_ylabel("\$x_1, y_1\$",labelpad = -20)
    ax2.set_ylabel("\$y_1\$")
    if i < M
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
    else
        ax1.set_xlabel("\$t\$"; labelpad = -20)
        ax2.set_xticks([-10, 0, 10])
        ax2.set_xticklabels(["-10", "", "10"])
        ax2.set_xlabel("\$x_1\$", labelpad = -20)
    end
end
add_identifiers!(fig, axs_for_labels)
fig.tight_layout(;pad=0.33)
wsave(plotsdir("9", "coupled_roessler"), fig)
