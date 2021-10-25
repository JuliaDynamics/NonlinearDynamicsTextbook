# 1d climate
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

using DynamicalSystems, PyPlot, Roots

αtan(T) = 0.5 - 0.2*tanh((T-263)/4)
dTdt(T, ε = 0.65, α=αtan, s= 1.0) = s*(1 - α(T)) - 1.6e-10 * ε * T^4
dTdt(T; ε = 0.65, α=αtan, s = 1.0) = dTdt(T, ε, α, s)

# Ts = 200:400.0
# plot(Ts, dTdt.(Ts))
# plot(Ts, dTdt.(Ts, 0.2))
# plot(Ts, dTdt.(Ts, 0.9))
# axhline(0)

fig = figure(figsize = (figx/2, figy))
Ts = 200:0.5:320.0
arrows = 210:10:300 |> collect
deleteat!(arrows, findfirst(isequal(260), arrows))
roots = Roots.find_zeros(dTdt, Ts[1], Ts[end])
plot(Ts, dTdt.(Ts), color = "C0", label = "\$dT/dt\$")
axhline(0; lw = 2.0, zorder = 1, color = "C2", ls = "--")
xlim(Ts[1], Ts[end])
ylim(-0.2, 0.2)
yticks([-0.1, 0, 0.1])
tight_layout()
for (i, r) in enumerate(roots)
    plot(r, 0, marker = "o", markeredgecolor = "k", markerfacecolor = iseven(i) ? "w" : "k",
    markersize = 15, mew = 2)
end
for r in arrows
    f = dTdt(r)
    x, dx = f > 0 ? (r - 5, 5) : (r+5, -5)
    ff = abs(1.2f)^2 + 0.01
    arrow(x, 0, dx, 0; color = "C1", width = ff, length_includes_head = false,
    head_width = 1.5ff, head_length = 100ff, zorder = 99
    )
end
legend()
xticks(200:40:320)
xlabel("\$T\$"; labelpad = -15)
tight_layout(pad=0.3)
wsave(plotsdir("2", "1dstatespace"), fig)
