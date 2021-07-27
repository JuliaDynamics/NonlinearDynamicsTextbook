using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

function coupled_roesslers(u, p, t)
    @inbounds begin
    a, b, c, ω1, ω2, k = p
    du1 = -ω1*u[2]-u[3]
    du2 = ω1*u[1] + a*u[2] + k*(u[5] - u[2])
    du3 = b + u[3]*(u[1] - c)
    du4 = -ω2*u[5]-u[6]
    du5 = ω2*u[4] + a*u[5] + k*(u[2] - u[5])
    du6 = b + u[6]*(u[4] - c)
    return SVector{6}(du1, du2, du3, du4, du5, du6)
    end
end

p0 = [0.15, 0.2, 10, 1-0.015, 1+0.015, 0.05]
u0 = rand(6)

croesslers = ContinuousDynamicalSystem(coupled_roesslers, u0, p0)

dt = 0.1
T = 1000
tr = trajectory(croesslers, T; Ttr = 100, dt)
t = 0:dt:T
L = length(t)÷10

fig = figure()
ax1 = fig.add_subplot(1,3,(1, 2))
ax2 = fig.add_subplot(1, 3, 3)
ax1.plot(t[1:L], tr[1:L, 1])
ax1.plot(t[1:L], tr[1:L, 4])
ax2.plot(tr[:, 1], tr[:, 4]; lw = 1.0)

# %% Lyapunovs
ks = range(0, 0.1; length = 200)

λs = zeros(length(ks), 4)

for (i, k) in enumerate(ks)
    @show i, k
    set_parameter!(croesslers, 6, k)
    λs[i, :] .= lyapunovspectrum(croesslers, 50000, 4; Ttr = 100)
end

fig, ax = subplots()
ax.plot(ks, λs)
ax.set_xlabel("\$k\$")
ax.set_ylabel("\$\\lambda_i\$")

fig.tight_layout(;pad = 0.25)
wsave(plotsdir("coupled_roess_lyap_chaotic"), fig)
