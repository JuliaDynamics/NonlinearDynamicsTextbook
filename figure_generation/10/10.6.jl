# Power grid/braess paradox
# From https://link.springer.com/article/10.1140/epjb/e2013-40469-4
# Nonlocal failures in complex supply networks by single link additions
# Dirk Witthaut & Marc Timme 
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random, OrdinaryDiffEq

using LinearAlgebra
P = 1.0
K = 1.03*P
α = 1.0*P
# Coefficient of power of the nodes, according to figure
power = Float64[-1, +1, +1, +1, -1, -1, +1, -1]
N = length(power) # number of nodes
# Weighted adjacency matrix, in units of K:
links = zeros(N, N)
links[1,2] = links[1,6] = links[1,7] = 1
links[2,1] = links[2,3] = 1
links[3,2] = links[3,4] = links[3,7] = 1
links[4,3] = links[4,5] = links[4,8] = 1
links[5,4] = links[5,6] = 1
links[6,5] = links[6,1] = links[6,8] = 1
links[7,1] = links[7,3] = 1
links[8,4] = links[8,6] = 1
@assert issymmetric(links)
# Colors of consumer/producer
PLUS = "C1"
MINUS = "C3"

function power_grid_f!(du, u, p, t)
    P, K, α, power, links = p
    N = size(u, 1)
    # φ is first column, ω is second column of `u`
    @inbounds for i in 1:N
        du[i,1] = u[i,2] # φdot = ω
        du[i,2] = power[i]*P - α*u[i,2] + 
            sum(K*links[i,j]*sin(u[j,1] - u[i,1]) for j in 1:N)
    end
end

u0 = zeros(N,2)
T = 20.0
saveat = 0:0.1:T

function rungrid!(ax, links)
    p = (P, K, α, power, links)
    pg = ODEProblem(power_grid_f!, u0, (0.0, T), p)
    sol = solve(pg, Tsit5(); saveat)
    # Phases timeseries
    φs = []
    for i in 1:N
        push!(φs, [u[i,1] for u in sol.u])
    end
    for i in 1:N
        c = power[i] > 0 ? PLUS : MINUS
        # ax.plot(sol.t, (φs[i] .- φs[1]), alpha = 0.9, color = c, lw = 2)
        ax.plot(sol.t, φs[i], alpha = 0.9, color = c, lw = 2)
    end
end

fig = figure()
rungrid!(gca(), links)


# %% plotting
# node coordinates
nx = [2.5, 2, 2.5, 0.5, 0, 0.5, 3, 1]
ny = [-1, 0, 1, 1, 0, -1, 0, 0]

fig, axs = subplots(2,3)

function plot_graph!(ax, nx, ny, power, links)
    colors = [p>0 ? "C1" : "C3" for p in power]
    colors = darken_color.(colors, 1.2)
    d = darken_color.(colors)
    ax.scatter(nx, ny; c = colors, s = 600, zorder = 98, edgecolors = d, linewidths = 2)
    for i in 1:N # loop over nodes
        t = power[i] > 0 ? "+" : "-"

        ax.annotate(t, (nx[i], ny[i]), fontsize = 24, ha = :center, va = :center, zorder = 99)
        # ax.annotate(i, (nx[i], ny[i]), fontsize = 24)
        # links
        for j in 1:N
            k = links[i,j]
            k == 0 && continue
            if abs(k) == 1
                color = BLACK
                lw = 3
            else
                color = "C0"
                lw = 9
            end
            ax.plot([nx[i],nx[j]], [ny[i],ny[j]]; color, lw, zorder = 1)
        end
    end
end

links2 = copy(links)
links2[3,4] = links2[4,3] = 2
links3 = copy(links)
links3[4,2] = links3[2,4] = 1
alllinks = (links, links2, links3)

for j in 1:3
    rungrid!(axs[2,j], alllinks[j])
    plot_graph!(axs[1,j], nx, ny, power, alllinks[j])
    axs[1,j].axis("off")
    axs[1,j].set_ylim(-1.4,1.4)
    axs[1,j].set_xlim(-0.4,3.4)
    y = 6
    axs[2,j].set_ylim(-y,y)
    axs[2,j].set_yticks([-y,y])
    if j > 1
        plt.setp(axs[2,j].get_yticklabels(), visible=false)
    end
    axs[2,j].set_xticks(0:6:18)
    axs[2,j].set_xlim(0,18)
    axs[2,j].set_xlabel("\$t\$",labelpad = -20)
    t = j == 1 ? "sync." : "no sync."
    axs[2,j].text(0.05, 0.8, t, transform = axs[2,j].transAxes)
end
axs[1,1].set_title("original")
axs[1,2].set_title("added capacity")
axs[1,3].set_title("added link")
axs[2,1].set_ylabel("\$\\phi_i\$", labelpad = -30)
add_identifiers!(fig, axs[2,:])
fig.tight_layout(pad=0.3)
wsave(plotsdir("10", "powergrid"), fig)