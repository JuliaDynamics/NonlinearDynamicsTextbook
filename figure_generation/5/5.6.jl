# %% Entropy comparison
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

he = Systems.henon()
N, q = 10000, 2
X = trajectory(he, N, Ttr = 10)
ε = 10 .^ (-5.0:0.1:1)
H = [genentropy(X, e; q) for e in ε]
C = correlationsum(X, ε)

fig, axs = subplots(2, 1; figsize = (figx/2, figy), sharex = true)

axs[1].plot(log.(ε), -H, COLORS[1])

x = -5:0
D = 1.22

axs[2].plot(log.(ε), log.(C), COLORS[1])

for (j, (ax, x)) in enumerate(zip(axs, (-5:0, -10:0)))
    Dl = D .* x
    Dl = Dl .- 2
    ax.axvspan(x[1], x[end], color = COLORS[4], alpha = 0.25)
    ax.grid("on")
    ax.plot(x, Dl, c = COLORS[3])
    ax.text(x[2] + 1.5 + 0.5*(j-1), Dl[2], "\$D_2\$ = $D", color = COLORS[3], size = 30)
end

xtick = -11:3:1
xtickl = string.(xtick)
xtickl[3] = ""
axs[2].set_xticks(xtick)
axs[2].set_xticklabels(xtickl)
axs[2].set_xlabel("\$\\log(\\varepsilon)\$"; labelpad = -20)
# axs[1].set_yticks()
axs[1].set_ylabel("\$-H_$q\$")
axs[2].set_ylabel("\$\\log(C)\$")
axs[1].set_yticks(0:-3:-10)
axs[2].set_yticks(0:-4:-12)
fig.tight_layout(pad=0.3)
wsave(plotsdir("5", "dimension"), fig)
