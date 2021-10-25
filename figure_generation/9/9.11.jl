# roessler chaotic phase diff
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random, OrdinaryDiffEq, LinearAlgebra

a = 0.2
b = 0.2
c = 5.7
amu = 0.02
ω1 = 1.0 + amu
ω2 = 1.0 - amu
k1 = 0.01
k2 = k1

y_init = [0.1, 0.2, 0., 0.11, 0.19, 0.1]

croesslers = Systems.coupled_roessler(y_init; ω1, ω2, a, b, c, k1, k2)
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
Δt = 0.01
T = 1000
t = 0:Δt:T
L = length(t)÷10

ks = [0.03, 0.04] # estimated from the Lyapunov exponents plot
bbox2 = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 0.8)

fig, axs = subplots(1,2; figsize = (0.66figx, figy))

function phase(x, y)
    f = atan.(y, x)
    for i in 2:length(f)
        if f[i] - f[i-1] < 0
            f[i:end] .+= 2π
        end
    end
    return f
end

for (i, k) in enumerate(ks)
    set_parameter!(croesslers, 6, k)
    set_parameter!(croesslers, 7, k)
    tr = trajectory(croesslers, T; Ttr = 200, Δt, diffeq...)
    x1 = tr[:, 1]
    x2 = tr[:, 2]
    y1 = tr[:, 4]
    y2 = tr[:, 5]
    f1 = phase(x1, x2)
    f2 = phase(y1, y2)
    axs[2].plot(t, f1 .- f2, label = "\$k=$(k)\$")

    if i == 1
        c = COLORS[1]
        axs[1].plot(x1[1:L], x2[1:L]; color = "0.5", linewidth = 2)
        j = 250
        axs[1].quiver(zeros(1), zeros(1), [x1[j+2]], [x2[j+2]], 
            color = c,
            scale = 1., scale_units = "xy", zorder = 9,
            width = 0.01
        )
        axs[1].plot([x1[j+1]], [x2[j+1]], marker = "o", zorder = 99, ms = 8,
        ls = "None", mec = c, mfc = :white,
            # scale = 1., scale_units = "xy", zorder = 9,
            # width = 0.02
        )
        θ = atan(x2[j], x1[j])
        r = norm([x1[j], x2[j]])
        a = matplotlib.patches.Wedge((0,0), 0.3r; 
        theta1 = 0, theta2 = 180θ/π, color = c, alpha = 0.5)
        axs[1].add_artist(a)
        axs[1].text(0.3x1[j-10], 0.3x2[j-10], "\$\\alpha\$"; color = c)
    end
end

axs[1].set_xlabel("\$x_1\$"; labelpad = -20)
axs[1].set_ylabel("\$x_2\$"; labelpad = -20)
axs[1].set_xticks([-10, 0, 10])
axs[1].set_xticklabels([-10, "", 10])
axs[2].set_xlabel("\$t\$"; labelpad = -20)
axs[2].set_ylabel("\$\\Delta \\alpha\$")
axs[2].set_xticks([0,t[end]])
axs[2].set_xlim(0,t[end])
axs[2].legend()
add_identifiers!(fig)
fig.tight_layout(;pad=0.33)
wsave(plotsdir("9", "coupled_roessler_phasediff"), fig)
