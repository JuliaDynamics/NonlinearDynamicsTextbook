using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using OrdinaryDiffEq, DynamicalSystems
using LinearAlgebra: mul!
import FFTW

# In the following `u` is the field in real space, while `y` is the Fourier
# transform of the field, which lives in inverse space.

b1 = 30 # total domain size
T = 500.0
dt = 0.25 # timestep for integration
saveat = 0:dt:T
dt = 0.25
dx = 0.2 # spatial discretization
xs = range(0, b1, step = dx) # space
u0 = @. cos(xs) + 0.1*sin(xs/8) + 0.01*cos((2π/b)*xs)
ks = Vector(FFTW.rfftfreq(length(u0))/dx) # conjugate space (wavenumbers)

forward_plan = FFTW.plan_rfft(u0)
y0 = forward_plan * u0
inverse_plan = FFTW.plan_irfft(y0, length(u0))
ik2 = -im .* ks ./ 2
k²_k⁴ = @. ks^2 - ks^4

ydummy = copy(y0)
udummy = copy(u0)
ksparams = (forward_plan, inverse_plan, udummy, ydummy, k²_k⁴, ik2)

function kse_spectral!(dy, y, p, t)
    forward_plan, inverse_plan, udummy, ydummy, k²_k⁴, ik2 = p
    y² = begin # nonlinear term
        mul!(udummy, inverse_plan, y) # create current u in x-space
        udummy .*= udummy # square current u
        mul!(ydummy, forward_plan, udummy) # transform to k-space
    end
    # KS equation in spectral space
    @. dy = y²*ik2 + k²_k⁴*y
    return nothing
end

prob = ODEProblem(kse_spectral!, y0, (0.0, T), ksparams)
@time sol = solve(prob, Tsit5(); saveat)

u = [inverse_plan*y for y in sol.u]
U = hcat(u...)

# %% Make figure
fig, axs = subplots(1,2)
pm = axs[1].pcolormesh(saveat, xs, U, cmap = "PRGn", vmin = -3, vmax = 3)
cb = colorbar(pm, ax = axs[1])
axs[1].set_xticks([0, saveat[end]])
axs[1].set_yticks([0, b1])
axs[1].set_xlabel("\$t\$"; labelpad = -20)
axs[1].set_ylabel("\$x\$"; labelpad = -20)

cb.set_ticks([-3, 3])
cb.set_label("\$u\$", labelpad = -40)

# load exponents
λs = []

N = 4000    # steps to run the QR normalization
Δt = 1.0    # timestep between normalizations
dx = 0.2    # spatial discretization
bs = 20:20:60

for b in bs
    try 
        N = 4000
        config = @strdict N dx Δt b
        file = wload(datadir("ksiva", savename("ksiva_spectrum", config, "jld2")))
        push!(λs, file["λs"])
    catch
        N = 1000
        config = @strdict N dx Δt b
        file = wload(datadir("ksiva", savename("ksiva_spectrum", config, "jld2")))
        push!(λs, file["λs"])
    end
end

N = 4000    # steps to run the QR normalization
for dx in (0.15, 0.1)
    b = 20
    config = @strdict N dx Δt b
    file = wload(datadir("ksiva", savename("ksiva_spectrum", config, "jld2")))
    push!(λs, file["λs"])
end

axs[2].axhline(0; color = "C3", lw = 1)

for i in 1:length(λs)
    l = λs[i]
    if i ≤ length(bs)
        axs[2].plot(1:length(l), cumsum(l); color = "C$(i-1)",
        label = "\$b=$(bs[i])\$")

        Δ = kaplanyorke_dim(sort(l; rev = true))
        axs[2].plot([Δ, Δ], [-25, 0]; color = "C$(i-1)", marker = "o", ms = 12,
        mec = "C$(i-1)", mfc = :white, lw = 1, zorder = 99)
    else
        axs[2].plot(1:length(l), cumsum(l); color = "C0",
        ls = i == length(bs) + 1 ? "--" : ":", zorder = 1)
    end
end

axs[2].set_ylim(-5, 5)
axs[2].set_xticks(0:40:120)
axs[2].set_xlim(0, 120)
axs[2].set_xlabel("\$n\$"; labelpad = -20)
axs[2].set_ylabel("\$\\Lambda_n\$", labelpad = -10)
axs[2].set_ylim(-5, 5)
axs[2].legend(ncol=1, loc = "upper right", fontsize = 24, handlelength = 1)
add_identifiers!(fig, axs)
fig.tight_layout(pad=0.3)
wsave(plotsdir("11", "ksiva"), fig)
