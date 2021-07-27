# %% Fractal dimension of chaotic attractors
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Statistics, DelimitedFiles

fig = figure()
ax_ts = fig.add_subplot(1,3,1)
ax_3d = fig.add_subplot(1,3,2; projection = "3d")
ax_fd = fig.add_subplot(1,3,3)

vostok = readdlm(projectdir("exercise_data", "17.csv"))
time = vostok[:, 1]
δT = vostok[:, 2]

# ax_ts.plot(time, δT; lw = 1.0, alpha = 0.5)

# Interpolate Vostok data
using Dierckx
# interpolation object:
spl = Spline1D(time, δT; k=3, bc="extrapolate")
# equi-spaced time vector:
t = range(minimum(time), maximum(time); length = 2000)

T = spl(t)

ax_ts.plot(t ./ 1000, T; lw = 1.0)
ax_ts.set_xlabel("time (kyr)")
ax_ts.set_ylabel("\$\\delta T\$ (K)")

ax_ts.set_xticks([0, 200, 400])
ax_ts.set_xticklabels([0, -200, -400])

# mi = selfmutualinfo(T, 1:1000)
# ax_fd.plot(1:1000, mi)
τ = 75 # inspect the above plot to see a good delay time

# plot attractor
R = embed(T, 3, τ)
ax_3d.plot(columns(R)...; lw = 1.0)
ax_3d.set_xticklabels([])
ax_3d.set_yticklabels([])
ax_3d.set_zticklabels([])
ax_3d.dist=8
ax_3d.text(5, 0, 7, "b"; bbox = bbox, zorder = 99, va = "top")

# estimate fractal dim
ds = 3:12
colors = matplotlib.cm.inferno(ds ./ (ds[end] + 2))
es = estimate_boxsizes(R; z = 0)

for (i, d) in enumerate(ds)
    A = embed(T, d, τ)
    Cs = boxed_correlationsum(A, es)
    j = findfirst(z -> z > 0, Cs)
    x, y = log.(es)[j:end], log.(Cs)[j:end]
    # x, y = log.(es), log.(Cs)
    # li, Δ = linear_region(x, y)
    c = colors[i, :]
    ax_fd.plot(x, y; color = c, lw = 1.5)
    # li = findall(z -> -4 ≤ z ≤ 0, x)
    _, Δ = linreg(x, y)
    @show (d, Δ)
    # ax_fd.plot(x[[1, end]], y[[1, end]],
    # zorder = 2, ms = 10, marker = "o", ls = "None", color = c, alpha = 0.75,
    # label = "d=$(d), Δ=$(round(Δ; digits=2))", )
end

# Show two regions
# ax_fd.axvspan(-3, -0.25; color = "C1", alpha = 0.1)
el1 = matplotlib.patches.Ellipse(
    (-1, -10), 3, 5; angle = -15, linewidth = 2, color = "C1",
    fill = false
)
ax_fd.add_patch(el1)
el2 = matplotlib.patches.Ellipse(
    (0.5, -5), 4, 6; angle = -20, linewidth = 2, color = "C2",
    fill = false
)
ax_fd.add_patch(el2)

ax_fd.set_xlabel("\$\\log(\\varepsilon)\$")
ax_fd.set_ylabel("\$\\log(C)\$";labelpad = -10)
ax_fd.set_yticks(0:-4:-12)
add_identifiers!(fig)
fig.tight_layout(;pad = 0.3)
fig.tight_layout(;pad = 0.33)
# wsave(plotsdir("vostok"), fig)
