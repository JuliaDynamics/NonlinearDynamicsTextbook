# %% mutual information
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

using DynamicalSystems, PyPlot, Random
Random.seed!(2586001)
lo, N = Systems.logistic(), 100
x = trajectory(lo, N-1)
y = trajectory(lo, N-1; Ttr = 1)
x .+= randn(N)/25; y .+= randn(N)/25

Hx = genentropy(Dataset(x), 0.1)
Hy = genentropy(Dataset(y), 0.1)
Hxy = genentropy(Dataset(x, y), 0.1)
m = Hx + Hy - Hxy # Eq. (7.1)
null = zeros(10000)
for i in 1:10000
  shuffle!(x); shuffle!(y);
  Hxy = genentropy(Dataset(x, y), 0.1)
  null[i] = Hx + Hy - Hxy
end

# plot stuff
fig, ax = subplots(;figsize = (figx/3, 1.25figy))
using Statistics
μ, σ = mean(null), std(null)
ax.hist(null, 50, label = "null pdf")
ax.axvline(μ, color = "C1", label = "\$\\mu\$", ls = "dashed")
ax.axvline(μ-3σ, color = "C2", label = "\$\\mu \\pm 3\\sigma\$", ls = "dashed")
ax.axvline(μ+3σ, color = "C2", ls = "dashed")
ax.axvline(m, color = "C3", label = "\$m\$")
ax.legend(
    bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
    ncol=2, mode="expand", borderaxespad=0, handlelength=1,
    fontsize = 26,
)

ax.set_yticks([])
ax.set_xlabel("mutual inform. (a. u.)")
fig.tight_layout(;pad = 0.25)
wsave(plotsdir("mutualinfo"), fig)
