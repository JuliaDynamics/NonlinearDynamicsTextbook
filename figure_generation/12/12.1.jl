# %% Sparse state space
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

using LinearAlgebra, Statistics

function sparse_f!(xnew, x, p, n)
    xnew[1] = 4*x[1]*(1-x[1])
    @inbounds for i in 2:length(x)
        xnew[i] = p[i]*x[i-1]
    end
    return
end

function sparse_jacob!(J, x, p, n)
    J[1,1] = 4 - 8x[1]
end

function make_sparse_system(D, p = [0.2*rand() + 0.87 for _ in 1:D-1])
    u0 = rand(D)
    # Optimized Jacobian in the form of tri-diagonal
    J = Tridiagonal(zeros(D-1), zeros(D), zeros(D-1))
    # Add entries of linear components, that will not change in time
    for i in 2:D; J[i, i-1] = p[i-1]; end
    # Adjust the first entry of the Jacobian, which depends on the state
    sparse_jacob!(J, u0, p, 0)
    sparse = DiscreteDynamicalSystem(sparse_f!, u0, p, sparse_jacob!, J)
end

sparse = make_sparse_system(400)
# A = trajectory(sparse, 1000; save_idxs = 1)
#
# plot(A[:, 1])

# Ds = 10:20:400
#
# λs = lyapunovspectrum(sparse, 10000, 10)
#
# Δs[i] = kaplanyorke_dim(sort!(λs; rev = true))

function sparse_lyapunov(config)
    D = config["D"]
    sparse = make_sparse_system(D)
    @time λs = DynamicalSystems.lyapunovspectrum(sparse, 10000, min(D, 20))
    Δ = kaplanyorke_dim(sort!(λs; rev = true))
    ret = Dict{String, Any}(config)
    @pack! ret = λs, Δ
    return ret
end

ret = sparse_lyapunov(Dict("D" => 10))

# data 1: kaplan yorke for increasing D
fig, axs = subplots(1,2; figsize = (figx, figy))

Ds = 10:10:200
Δs = zeros(length(Ds))

for (i, D) in enumerate(Ds)
    @show D
    file, _ = produce_or_load(datadir("logisticdrift"), Dict("D" => D), sparse_lyapunov; 
    tag = false)
    # sparse = make_sparse_system(D)
    # @time λs = DynamicalSystems.lyapunovspectrum(sparse, 10000, 10)
    # Δs[i] = kaplanyorke_dim(sort!(λs; rev = true))
    Δs[i] = file["Δ"]
end
# plot 1


axs[1].plot(Ds, Δs)
axs[1].set_xlabel("\$D\$"; labelpad = -10)
axs[1].set_ylabel("\$\\Delta^{(L)}\$"; labelpad = -5)
axs[1].set_yticks(1:4)
axs[1].set_xticks(10:40:200)


# %% data 2: perturbation growth lala
# fig, axs = subplots(1,2; figsize = (figx, figy))
using Random
# Random.seed!(77670011)
Ds = 10 .^ (1:4)
for (d, D) in enumerate(Ds)
    @show D
    sparse = make_sparse_system(D)


    S = 1000 # samples
    n = 50 # max iteration
    lD = zeros(n, S)

    for j in 1:S
        Q0 = normalize!(rand(D, 1))
        while any(isnan, Q0)
            Q0 = normalize!(rand(D, 1))
        end
        Q0 .* 1e-6
        tinteg = tangent_integrator(sparse, Q0; u0 = rand(D))
        g = zeros(n)
        for i in 1:n
            step!(tinteg)
            g[i] = norm(get_deviations(tinteg))
        end

        # create instantaneous exponential growth rates
        lD[:, j] = log.(g) ./ (1:n)
    end

    lDμ = vec(mean(lD; dims = 2))
    lDσ = vec(std(lD; dims = 2))
    axs[2].plot(1:n, lDμ; color = "C$((0,1,3,4)[d])", label = "\$D=$D\$")
    # axs[2].fill_between(1:n, lDμ .- lDσ, lDμ .+ lDσ; color = "C$(d)", alpha = 0.5)
end

# add lyapunov exponent indication
axs[2].axhline(log(2); ls = "dashed", color = "C2", lw = 1.5)
axs[2].set_xlabel("\$n\$"; labelpad = -10)
axs[2].set_xticks(0:16:48)
axs[2].set_ylabel("\$\\lambda_\\mathrm{local}\$")
axs[2].legend()
axs[2].set_ylim(-0.1, 0.75)
add_identifiers!(fig)
fig.tight_layout(;pad = 0.3)
# wsave(plotsdir("sparse_space"), fig)
