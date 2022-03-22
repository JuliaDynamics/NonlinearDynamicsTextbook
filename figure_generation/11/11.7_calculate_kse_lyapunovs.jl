# This code calculates the lyapunov spectrum of the KS equation 
# by transforming it in spectral space.
# It uses very efficient handling of Fourier transforms.
# It saves exponents to be used by the figure generating scripts.
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using OrdinaryDiffEq
import FFTW, LinearAlgebra

# In the following `u` is the field in real space, while `y` is the Fourier
# transform of the field, which lives in inverse space.
# `w` is a matrix that has the first column the state `y` and the remaining the
# deviation vectors `δy`.

function kse_spectral_and_tangent!(dw, w, p, t)
    @unpack  ks, forward_plan, inverse_plan, yd, ud, ud2, ud3, yd2, ik2, k²_k⁴ = p
    y = view(w, :, 1)
    dy = view(dw, :, 1)

    # KS equation in spectral space
    LinearAlgebra.mul!(ud, inverse_plan, y) # create current u in x-space
    @. ud3 = ud*ud
    y² = LinearAlgebra.mul!(yd, forward_plan, ud3) # transform to k-space
    @. dy = - y²*ik2 + k²_k⁴*y

    # Equations for tangent space in spectral space (loop over each deviation vector)
    for i in 2:size(w, 2)
        δy = view(w, :, i)
        dδy = view(dw, :, i)
        ud2 .= ud # keeps track for tangent
        LinearAlgebra.mul!(ud2, inverse_plan, δy)
        ud2 .*= ud
        spectral_tan = LinearAlgebra.mul!(yd2, forward_plan, ud2)
        @. dδy = k²_k⁴*δy + im*ks*spectral_tan
    end
    return nothing
end

function orthonorm_spectral(D, M = D) # Orthonormal vectors converted to spectral space
    M > D && throw(ArgumentError("M must be ≤ D"))
    Δ = Matrix(LinearAlgebra.qr(rand(D, M)).Q)
    FFTW.rfft(Δ)
end

function tannorm(u::AbstractMatrix, t)
    s = size(u)[1]
    x = 0.0
    for i in 1:s; @inbounds x += abs2(u[i, 1]); end
    return sqrt(x/length(x))
end
tannorm(u, t) = abs(u)

# %% Lyapunov code
# This function uses some convenience methods from DynamicalSystems.jl
# based on the `tangent_integrator` format (state is a matrix. first column
# is actual system state, all other columns are deviation vectors)
import ProgressMeter
function kse_lyapunovs_spectral(integ, N, Δt::Real)
    progress = ProgressMeter.Progress(N; desc = "KSE Lyapunov Spectrum: ", dt = 1.0)
    M = size(integ.u, 2) - 1 # number of Ls
    λ = zeros(M)
    forward_plan = integ.p.forward_plan
    inverse_plan = integ.p.inverse_plan
    W = zeros(length(inverse_plan*integ.u[:, 1]), M) # for use in buffer
    t0 = integ.t
    ud = zeros(size(W, 1))
    yd = integ.u[:, 1]

    for n in 1:N
        step!(integ, Δt)
        # Get deviations in real space 
        for i in 1:M
            LinearAlgebra.mul!(ud, inverse_plan, view(integ.u, :, i+1))
            W[:, i] .= ud
        end
        # Perform a (buffered) QR decomposition
        Q, R = LinearAlgebra.qr!(W)
        # Keep track of LEs
        for j in 1:M
            @inbounds λ[j] += log(abs(R[j,j]))
        end
        # Set the new deviations (also convert back to spectral)
        for i in 1:M
            ud .= Q[:, i]
            LinearAlgebra.mul!(yd, forward_plan, ud)
            integ.u[:, i+1] .= yd
        end
        u_modified!(integ, true) # Ensure that DiffEq knows we did something!
        ProgressMeter.update!(progress, n)
    end
    λ ./= (integ.t - t0)
    return λ
end

# %% Actually compute the exponents
function produce_kse_lyaps(config)
    @unpack b, dx, N, Δt, kport = config
    xs = range(0, b; step = dx) # space
    u0 = @. cos(xs) + 0.1*sin(xs/8) + 0.01*cos((2π/b)*xs)
    ks = Vector(FFTW.rfftfreq(length(u0))/dx) # conjugate space (wavenumbers)

    forward_plan = FFTW.plan_rfft(u0)
    y0 = forward_plan * u0
    inverse_plan = FFTW.plan_irfft(y0, length(u0))
    ik2 = -im .* ks ./ 2
    k²_k⁴ = @. ks^2 - ks^4

    ud = copy(u0)
    ud2 = copy(u0)
    ud3 = copy(u0)
    yd = copy(y0)
    yd2 = copy(y0)
    ksparams = @ntuple forward_plan inverse_plan ks ud ud2 ud3 yd yd2 k²_k⁴ ik2

    D = length(u0)
    M = round(Int, kport*D)
    δy0 = orthonorm_spectral(D, M)
    w0 = hcat(y0, δy0) 
    prob = ODEProblem(kse_spectral_and_tangent!, w0, (0.0, 100.0), ksparams)
    integ = init(prob, Tsit5(); save_everystep = false, internalnorm = tannorm)

    @time λs = kse_lyapunovs_spectral(integ, N, Δt)
    return @strdict λs b dx N Δt kport
end

b = 20      # length of spatial domain
N = 4000    # steps to run the QR normalization
Δt = 1.0    # timestep between normalizations
dx = 0.2    # spatial discretization

for b in 20:20:60
    if b == 20
        kport = 0.5
    else
        kport = 0.3
    end
    @show b
    config = @strdict N dx Δt b kport
    file, path = produce_or_load(
        datadir("ksiva"), config, produce_kse_lyaps; 
        prefix = "ksiva_spectrum", ignores = ["kport"], storepatch = false,
    )
end

b = 20

for dx in (0.2, 0.15, 0.1)
    kport = 0.5
    @show dx
    config = @strdict N dx Δt b kport
    file, path = produce_or_load(
        datadir("ksiva"), config, produce_kse_lyaps; 
        prefix = "ksiva_spectrum", ignores = ["kport"], storepatch = false,
    )
end