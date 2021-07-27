using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

# %% The very first code snippet of the course
using DynamicalSystems # load the library

function lorenz_rule(u, p, t)
    σ, ρ, β = p
    x, y, z = u
    dx = σ*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z
    return SVector(dx, dy, dz) # Static Vector
end

p  = [10.0, 28.0, 8/3] # parameters: σ, ρ, β
u₀ = [0, 10.0, 0]      # initial state
# create an instance of a `DynamicalSystem`
lorenz = ContinuousDynamicalSystem(lorenz_rule, u₀, p)

T  = 100.0 # total time
dt = 0.01  # sampling time
A  = trajectory(lorenz, T; dt)

# %% Animate evolution of trajectories in the Standard map
using InteractiveDynamics
using DynamicalSystems
using GLMakie

# Standard map trajectories
ds = Systems.standardmap(; k = 1.0)
u0s = [[θ, p] for θ ∈ 0:2π for p ∈ 0:2π]
lims = ((0, 2π), (0, 2π))

figure, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 1000, lims,
)
main.xlabel = "θ"
main.ylabel = "p"

# %% Lorenz system trajectories
using OrdinaryDiffEq: Tsit5, Vern9
ds = Systems.lorenz()

u0s =  [[10,20,40.0] .+ rand(3) for _ in 1:6]

diffeq = (alg = Tsit5(), dtmax = 0.01)

figure, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 1000, diffeq, colors = to_color.(COLORS)
)

# %% Add poincare plane
import AbstractPlotting
using AbstractPlotting.MakieLayout: LAxis
o = AbstractPlotting.Point3f0(-25, 0, 0)
w = AbstractPlotting.Point3f0(50, 0, 50)
p = AbstractPlotting.FRect3D(o, w)

# These lines are necessary for transparent planes
a = AbstractPlotting.RGBAf0(0,0,0,0)
c = AbstractPlotting.RGBAf0(0.2, 0.2, 0.8, 1.0)
img = AbstractPlotting.ImagePattern([c a; a c]);
AbstractPlotting.mesh!(main, p; color = img);

# %% Plot Poincare sos
psosplot = layout[:, 2] = LAxis(figure)
psos = poincaresos(ds, (2, 0.0), 2000.0)
AbstractPlotting.scatter!(psosplot, psos[:, 1], psos[:, 3])

display(figure)

# %% Henon heiles system trajectories
ds = Systems.henonheiles()

u0s = [[0.0, -0.25, 0.42081, 0.0],
[0.0, 0.1, 0.5, 0.0],
[0.0, -0.31596, 0.354461, 0.0591255]]

diffeq = (alg = Vern9(), dtmax = 0.01)
idxs = (1, 2)

figure, main, layout, obs = interactive_evolution(
    ds, u0s; idxs, tail = 2000, diffeq, colors = COLORS,
)
main.figure[AbstractPlotting.Axis][:names, :axisnames] = ("q₁", "q₂", "p₂")

# %% Poincare brainscanning application
using GLMakie, DynamicalSystems, InteractiveDynamics
using OrdinaryDiffEq

ds = Systems.henonheiles()
diffeq = (alg = Vern9(),)
u0s = [
    [0.0, -0.25, 0.42081, 0.0],
    [0.0, 0.1, 0.5, 0.0],
    [0.0, -0.31596, 0.354461, 0.0591255]
]
trs = [trajectory(ds, 10000, u0; diffeq...)[:, SVector(1,2,3)] for u0 ∈ u0s]
for i in 2:length(u0s)
    append!(trs[1], trs[i])
end

# Inputs:
j = 1 # the dimension of the plane
tr = trs[1]

brainscan_poincaresos(tr, j)
