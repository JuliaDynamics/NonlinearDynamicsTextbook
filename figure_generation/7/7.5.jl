# %% Surrogate application (From Lancaster)
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

using TimeseriesSurrogates, Statistics, Random

ro = Systems.roessler(ones(3); a = 0.165, b = 0.2, c = 10.0)

tr = trajectory(ro, 500; Δt = 0.1, Ttr = 500)
x = tr[:, 1]
x ./= std(x)
x .+= randn(length(x))*std(x)*0.1

# generate autoregressive process
Random.seed!(77163)
η = randn(5000)
s = ones(5000)
for n in 4:5000
    s[n] = 1.625s[n-1] - 0.284s[n-2] - 0.355s[n-3] + η[n] - 0.96η[n-1]
end
s ./= std(s)

# Do the surrogate calculation
εro, εma = std.((x, s))./4
algs = [RandomFourier(), AAFT()]
names = ["FT", "AAFT"]
sgx = [surrogenerator(x, m) for m in algs]
sgs = [surrogenerator(s, m) for m in algs]
τx = estimate_delay(x, "ac_zero")
τs = estimate_delay(s, "ac_zero")
Cx = grassberger(embed(x, 3, τx))
Cs = grassberger(embed(s, 4, τs))
A = length(algs)

xboxes = []
sboxes = []
for i in 1:A
    sx, ss = sgx[i], sgs[i]
    bx, bs = [], []
    for j in 1:100
        X = embed(sx(), 3, τx)
        S = embed(ss(), 4, τs)
        push!(bx, grassberger(X))
        push!(bs, grassberger(S))
    end
    push!(xboxes, bx)
    push!(sboxes, bs)
end

# %% Plot everything
fig = figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,4,3)
ax3 = fig.add_subplot(1,4,4)
axs = [ax1, ax2, ax3]

axs[1].plot(x .+ 2,  label = "Rössler", lw = 1.5)
axs[1].plot(s .- 2,label = "ARMA", lw = 1.5, color = "C2")
axs[1].set_xlim(0, 1000)
axs[1].set_ylim(-5, 5)
axs[1].set_yticks(-4:2:4)
leg = axs[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
           ncol=2, mode="expand", borderaxespad=0., fontsize = 26)
for line in leg.get_lines()
    line.set_linewidth(4.0)
end

for (j, boxes) in enumerate((xboxes, sboxes))
    c = j == 1 ? "C0" : "C2"
    ax = axs[j+1]
    for (i, b) in enumerate(boxes)
        ax.boxplot([b]; positions = [i],
        patch_artist=true, boxprops=Dict(:facecolor=>c, :color=>c),
        medianprops=Dict(:color=>"w"), flierprops=Dict(:markeredgecolor=>c))
    end
end

l = [1-0.2, A+0.2]
axs[2].plot(l, fill(Cx, 2), color = "C0", ls = "dashed")
# axs[2].plot(l, fill(2.75, 2), color = "C0", ls = "dashed")
axs[3].plot(l, fill(Cs, 2), color = "C2", ls = "dashed")
# axs[2].set_ylim(2.7, 3)
axs[2].set_yticks(2.6:0.2:3.0)
axs[2].set_title("\$\\Delta^{(C)}\$, Rössler")
axs[3].set_title("\$\\Delta^{(C)}\$, ARMA")
axs[3].set_yticks(3.4:0.3:4.1)

for ax in axs[2:3]
    ax.grid(false; axis = "x")
    ax.set_xticks(1:A)
    ax.set_xticklabels(names, rotation = 0, size = 20)
end

fig.tight_layout(pad=0.3)
fig.subplots_adjust(top = 0.85, left = 0.05, bottom = 0.1, right = 0.98, wspace = 0.3)
wsave(plotsdir("7", "surrogates"), fig)