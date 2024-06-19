# orbit evolution on a 2D torus
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("theme.jl"))
using DynamicalSystems
import GLMakie
using OrdinaryDiffEq

function quasiperiodic_f(u, p, t)
    # here we make the frequency ratio a state variable, because the
    # `interactive_trajectory` application doesn't allow different parameters for different
    # initial conditions
    ω = u[3]
    θdot = ω
    φdot = 1.0
    return SVector(θdot, φdot, 0.0)
end

frequencies² = [7, 9]
u0s = [[0, 0, √x] for x ∈ frequencies²]
diffeq = (alg = Tsit5(), adaptive = false, dt = 0.01)
ds = ContinuousDynamicalSystem(quasiperiodic_f, [0.0, 0.0, 0.0], nothing; diffeq)

# Here we define a projection from frequencies into a torus in 3D
R = 2.0
r = 1.0

function torus(u)
    θ, φ = u
    x = (R + r*cos(θ))*cos(φ)
    y = (R + r*cos(θ))*sin(φ)
    z = r*sin(θ)
    return SVector(x, y, z)
end
torusx(u) = ((θ, φ) = u; (R + r*cos(θ))*cos(φ))
torusy(u) = ((θ, φ) = u; (R + r*cos(θ))*sin(φ))
torusz(u) = ((θ, φ) = u; r*sin(θ))

observables = [torusx, torusy, torusz] # inefficient but oh well

lims = ((-R-r, R+r), (-R-r, R+r), (-3r, 3r))
plotkwargs = [(linewidth = i^2*2.0,) for i in 1:length(u0s)]

fig, obs = interactive_trajectory(
    ds, u0s; tail = 20000, colors = COLORSCHEME[1:length(u0s)],
    lims, plotkwargs, idxs = observables
)

main = GLMakie.content(fig[1,1][1,1])

freqs = collect("√"*string(s) for s in frequencies²)
main.title = "Frequency ratios: "*join(freqs, ", ")

# Create a wireframe for the torus
U = range(-π, π; length = 200)
V = range(-π, π; length = 50)

x1 = [(R + r*cos(θ))*cos(φ)      for θ in U, φ in V]
y1 = [(R + r*cos(θ))*sin(φ)      for θ in U, φ in V]
z1 = [r*sin(θ)                   for θ in U, φ in V]
GLMakie.wireframe!(main, x1,y1,z1;
    color = :black, linewidth = 0.01, transparency = true,
)

display(fig)

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)

# %% Clean figure for torus in book
GLMakie.hidedecorations!(main)
GLMakie.hidespines!(main)