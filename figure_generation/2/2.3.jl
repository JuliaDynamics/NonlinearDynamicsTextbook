# Brusselator figure
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

function brusselator_f(x, p, t)
    u, v = x
    a, b = p
    du = 1 - (b+1)*u + a*u^2*v
    dv = b*u - a*u^2*v
    return SVector(du, dv)
end

a = 0.3; b = 1.5
p0 = [a, b]
u0 = [2.0, 2.0]

ds = ContinuousDynamicalSystem(brusselator_f, u0, p0)

Δt = 0.001
T = 100
tr = trajectory(ds, T; Δt, Ttr = 100)
t = 0:Δt:T
u, v = columns(tr)

period = estimate_period(u, :zerocrossing, t)
ip = round(Int, period/Δt)

fig = figure(figsize = (figx/2, figy))
ax = gca()

sc = ax.scatter(u[1:ip], v[1:ip]; c = range(0, 2π, length = ip), cmap = :twilight, s = 20)
cb = colorbar(sc)
cb.set_ticks([0, 2π])
cb.set_ticklabels(["\$0\$", "\$2\\pi\$"])
cb.set_label("\$\\phi\$", labelpad = -30)
ax.set_xticks([0.5, 2.0])
ax.set_yticks([3, 5])
ax.set_xlabel("\$u\$", labelpad = -20)
ax.set_ylabel("\$w\$", labelpad = -10)
ax.plot(1, b/a; marker = "o", mec = "C2", mew = 2,
markersize = 10, mfc = :white)

r = ip ÷ 10
ax.scatter(u[1:r:ip], v[1:r:ip]; color = "C2", marker = "D", s = 150, edgecolors="white")

fig.tight_layout(pad=0.3)

# wsave(plotsdir("2", "brusselator"), fig)