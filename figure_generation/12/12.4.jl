# %% Basin stability of magnetic pendulum
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot
using LinearAlgebra, Statistics

d, α, ω = 0.3, 0.2, 0.5
xg = yg = range(-3, 3, length = 300)
ma = Systems.magnetic_pendulum(;γs=[1,1,0.2], d, α, ω)

# This function is faster than basins, because it only picks some random points
function random_statespace_fraction(∫, xg, yg, T = 10000)
    c = 0 # count of states converged to attractor of interest
    for i in 1:T
        x = rand()*(maximum(xg)-minimum(xg)) + minimum(xg)
        y = rand()*(maximum(yg)-minimum(yg)) + minimum(yg)
        reinit!(∫, SVector(x, y, 0, 0))
        step!(∫, 100.0)
        while ∫.u[3]^2 + ∫.u[4]^2 > 1e-3
            step!(∫)
        end
        s = SVector(∫.u[1], ∫.u[2])
        dmin, k = findmin([(s-m)⋅(s-m) for m in ma.f.magnets])
        if k == 3; c += 1; end
    end
    return c/T
end

γs = 0:0.01:1
Fs = zero(γs)
Fσs = zero(γs)
λs = zero(γs)
∫ = integrator(ma)

# TODO: Use `produce_or_load` here.
for (i, γ) in enumerate(γs)
    @show (i, γ)
    ∫.p.γs[3] = γ
    allfracs = [random_statespace_fraction(∫, xg, yg, 1000) for i in 1:5]
    Fs[i] = mean(allfracs)
    Fσs[i] = std(allfracs)
    # TODO: This is slightly inaccurate; the fixed point is not exactly
    # on the magnetic. But it shouldn't have much of an impact for
    # final figure.
    J = Array(ma.jacobian(SVector(ma.f.magnets[3]..., 0, 0), ma.p, 0))
    eee = eigvals(J)
    λs[i] = maximum(real.(eee))
end

savedict = @strdict λs γs Fs
wsave(datadir("magnetic_basin_stability.jld2"), savedict)


# %%
data = wload(datadir("magnetic_basin_stability.jld2"))
@unpack Fs, λs, γs = data
fig, axs = subplots(1, 3)
axs[1].plot(γs, Fs; label = "\$F\$")
axs[1].fill_between(γs, Fs .- Fσs, Fs .+ Fσs; color = "C0", alpha = 0.5)
# axs[1].errorbar(γs, Fs; label = "\$F\$", yerr = Fσs)
# axs[1].plot(γs, Fs; label = "\$F\$")
axs[1].plot(γs, λs; label = "\$\\lambda_1\$")
axs[1].legend()
axs[1].set_xlabel("\$\\gamma_3\$"; labelpad = -20)
axs[1].set_xticks([0, 1])
axs[1].set_xlim(0,1)
axs[1].set_yticks(-0.1:0.1:0.3)

# perform 2 nice detailed plots
LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])
γplots = (0.7, 0.25)
basins_pre = att_pre = nothing
for (i, γ) in enumerate(γplots)
    ma.p.γs[3] = γ
    basins, attractors = basins_of_attraction((xg, yg), ma)
    if i == 2
        match_attractors!(basins_pre, att_pre, basins, attractors)
    end
    axs[i+1].pcolormesh(xg, xg, basins'; cmap, shading = "gouraud")
    axs[i+1].set_xlabel("\$x\$")
    axs[i+1].set_ylabel("\$y\$")
    axs[i+1].set_xticks([])
    axs[i+1].set_yticks([])
    basins_pre = basins
    att_pre = attractors
end

add_identifiers!(fig)
fig.tight_layout(;pad = 0.33)
wsave(plotsdir("12", "basin_stability"), fig)
