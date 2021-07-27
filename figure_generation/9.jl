# %%
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot
localmaxima(x) = findall(i -> x[i-1] < x[i] && x[i] > x[i+1], 2:length(x)-1) .+ 1

duffing = Systems.duffing(; ω = 2.0)

as = 1:0.1:100
ω = duffing.p[1]
T = 500.0
fig = figure()
ax = gca()
u0 = trajectory(duffing, 1000)[end]

for (i, as_) in enumerate((as, reverse(as)))
    for a in as_
        set_parameter!(duffing, 2, a)
        tr = trajectory(duffing, T, u0; Ttr = 100)
        global u0 = tr[end]
        x = tr[:, 1]
        j = localmaxima(x)
        n = length(j)
        ax.plot(fill(a, n), x[j]; ls = "None", ms = 1,
        color = "C$(i-1)", marker = "o", alpha = 0.5)
    end
end
ax.set_xlabel("driving amplitude \$a\$")
ax.set_ylabel("local maxima")
fig.tight_layout(;pad = 0.25)
