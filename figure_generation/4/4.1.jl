# Bifurcation diagram for 1D climate model
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

include(srcdir("simple_continuation.jl"))

using ForwardDiff, PyPlot
αtan(T) = 0.5 - 0.2*tanh((T-263)/4)
dTdt(T, ε = 0.65, α=αtan, s= 1.0) = s*(1 - α(T)) - 1.6e-10 * ε * T^4
dTdt(T; ε = 0.65, α=αtan, s = 1.0) = dTdt(T, ε, α, s)
d²Tdt²(T, ε) = ForwardDiff.derivative(T -> dTdt(T, ε), T)
climate_f(u, p, t) = SVector(dTdt(u[1], p))
climate_J(u, p, t) = SMatrix{1,1}(d²Tdt²(u[1], p))

u0 = SVector(225.0)
p0 = 0.9
ds = ContinuousDynamicalSystem(climate_f, u0, p0, climate_J)

f, J = bifurcation_rule_form(ds)

pmin = 0.3
pmax = 1

dp0 = 0.05
dx0 = SVector(0.1)

xs, ps, stability = continuation(f, J, u0, p0;
    pmin, pmax, dp0, dx0, N = 1500
)

xs = [x[1] for x in xs]


# %% make figure
stable1 = findfirst(!, stability)
unstable1 = findlast(!, stability)

fig = figure(figsize = (figx/2, figy))
ax = gca()

plot(ps[1:stable1], xs[1:stable1]; c = "C2", ls = "solid", label="stable")
plot(ps[stable1+1:unstable1], xs[stable1+1:unstable1]; c = "C2", ls = "dashed", label="unstable")
plot(ps[unstable1+1:end], xs[unstable1+1:end]; c = "C2", ls = "solid")

# scatter bifurcation points
ax.plot([ps[stable1], ps[unstable1]], [xs[stable1], xs[unstable1]], ls = "None",
		marker="o", ms = 12, color = "k", mew = 1, mec = "k",
        fillstyle = "top", markerfacecoloralt= "w"
)

# plot arrows
_a = 300
arrow1ε = ps[_a:stable1] .- 0.03
arrow1T = xs[_a:stable1] .- 5
plot(arrow1ε, arrow1T; color = "C1", lw = 4.0)
arrow(arrow1ε[end], arrow1T[end], 0, 50; 
width = 0.004, head_length = 5.0, color = "C1", zorder = 99)

_b = 200
arrow2ε = ps[unstable1:unstable1+_b] .+ 0.03
arrow2T = xs[unstable1:unstable1+_b] .+ 5

plot(arrow2ε, arrow2T; color = "C0", lw = 4.0)
arrow(arrow2ε[1], arrow2T[1], 0, -50; 
width = 0.004, head_length = 5.0, color = "C0", zorder = 99)

# scatter(ps[stability], xs[stability]; c = "C3", legend = "stable")
# scatter(ps[stability], xs[stability]; c = "C3", legend = "stable")
legend()
xticks(0.3:0.2:0.9)
yticks(200:50:350)
ylabel("\$T^*\$", labelpad = -20)
xlabel("\$\\varepsilon\$"; labelpad=-20)
# yticks(230:50:310)
# ylim(210, 320)
# ax.xaxis.set_label_coords(0.4, -0.025)
fig.tight_layout(pad=0.3)
wsave(plotsdir("4", "bif_example"), fig)