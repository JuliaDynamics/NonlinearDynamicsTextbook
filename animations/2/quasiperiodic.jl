# orbit evolution on a 2D torus
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveDynamics

using OrdinaryDiffEq

R = 2.0
r = 1.0

function torus(u)
    θ, φ = u
    x = (R + r*cos(θ))*cos(φ)
    y = (R + r*cos(θ))*sin(φ)
    z = r*sin(θ)
    return SVector(x, y, z)
end

function quasiperiodic_f(u, p, t)
    # here we make the frequency ratio a state variable, because the
    # interactive_evolution application doesn't allow different parameters for different
    # initial conditions
    ω = u[3]
    θdot = ω
    φdot = 1.0
    return SVector(θdot, φdot, 0.0)
end

ds = ContinuousDynamicalSystem(quasiperiodic_f, [0.0, 0.0, 0.0], nothing)
frequencies² = [7, 9]
u0s = [[0, 0, √x] for x ∈ frequencies²]
diffeq = (alg = Tsit5(), dtmax = 0.01)

lims = ((-R-r, R+r), (-R-r, R+r), (-3r, 3r))

plotkwargs = [(linewidth = i^2*1.0,) for i in 1:length(u0s)]

fig, obs = interactive_evolution(
    ds, u0s; tail = 20000, diffeq, colors = COLORSCHEME[1:length(u0s)], 
    transform = torus, lims, plotkwargs
)

main = content(fig[1,1])

freqs = collect("√"*string(s) for s in frequencies²)
main.title = "Frequency ratios: "*join(freqs, ", ")

# Create a wireframe for the torus
# angles = [SVector(θ, φ) for θ ∈ range(0, 2π; length = 100) for φ ∈ range(0, 2π; length = 100)]
# coords = torus.(angles)
# wireframe!(main, coords)


U = range(-π, π; length = 200)
V = range(-π, π; length = 50)

x1 = [(R + r*cos(θ))*cos(φ)      for θ in U, φ in V]
y1 = [(R + r*cos(θ))*sin(φ)      for θ in U, φ in V]
z1 = [r*sin(θ)                   for θ in U, φ in V]
wireframe!(main, x1,y1,z1, shading = false, color = COLORS[3], linewidth = 0.1)


# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
