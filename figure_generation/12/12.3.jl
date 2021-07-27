# %% Rate dependent tipping illustration
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

# For simplicity, this code continues the code of 4/4.1.jl 
# and loads from there the `stable1, unstable1`, `ps` , `xs`, etc.

fig = figure(figsize = (figx/2, figy))
axbif = gca()

axbif.plot(ps[1:stable1], xs[1:stable1]; c = "C2", ls = "solid", label="stable")
axbif.plot(ps[stable1+1:unstable1], xs[stable1+1:unstable1]; c = "C2", ls = "dashed", label="unstable")
axbif.plot(ps[unstable1+1:end], xs[unstable1+1:end]; c = "C2", ls = "solid")

axe = axbif.inset_axes(;bounds = [0.6, 0.6, 0.35, 0.35], transform = axbif.transAxes)
axbif.set_ylim(210, 370)

function rate_dependent_system(du, u, p, t)
    Δε, Δt, t0 = p
    T, ε = u
    dεdt = (-(t - t0)*Δε / Δt^2) * exp(-(t-t0)^2/(2*Δt^2))
    du[1] = dTdt(T, ε)
    du[2] = dεdt
    nothing
end

t0 = 2500.0

for (i, p) in enumerate([(0.1, 40.0), (0.15, 40.0), (0.1, 60.0)])
    p0 = [p..., t0] # if Δε is 0, it doesn't matter
    u0 = [274, 0.78] # lol gotta rework this because ε goes beyond 1.

    ds = ContinuousDynamicalSystem(rate_dependent_system, u0, p0)

    maxt = 5000
    tr = trajectory(ds, maxt; dt = 1.0)
    T, ε = columns(tr)
    c = ("C0", "C1", "C3")[i]
    axbif.plot(ε, T; color = c)
    axe.plot(0:maxt, ε; color = c)
end
axe.set_xlim(2300, 2700)
axe.set_xlabel("\$t\$")
axe.set_xticks([])
axe.set_ylabel("\$\\epsilon\$"; labelpad = -20)
axe.grid(false)
axbif.set_xticks(0.3:0.2:0.9)
axbif.set_yticks(210:50:370)
axbif.set_xlabel("\$\\epsilon\$"; labelpad = -20)
axbif.set_ylabel("\$T\$"; labelpad = -20)
fig.tight_layout(pad = 0.25)
# wsave(plotsdir("rate_dependent"), fig)
