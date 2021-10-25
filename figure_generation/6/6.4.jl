# %% Using Cao's method to estimate embedding
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random


lo = Systems.lorenz([0, 10, 0.0])

tr = trajectory(lo, 1000; Ttr=10)
x, y, z = columns(tr)
fig = figure(figsize=(figx/2, figy));
w = x
ψ = z .- y

for (i, (s, l)) in enumerate(zip((w, ψ), ("\$w=x\$", "\$w=z-y\$")))
    τ = estimate_delay(s, "mi_min")
    Ed = delay_afnn(s, τ, 2:7)
    plot(2:7, Ed, marker = "o", label = l, color = "C$([0, 2][i])")
end

xlabel("\$d\$"; labelpad = -10)
ylabel("\$E_{d}\$")
legend(title="measurement")
fig.tight_layout(pad = 0.4)
fig.subplots_adjust(left = 0.15, bottom = 0.15)
wsave(plotsdir("6", "caodemonstration"), fig)