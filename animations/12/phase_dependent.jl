using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using GLMakie, InteractiveDynamics, DynamicalSystems, OrdinaryDiffEq, LinearAlgebra
diffeq = (alg = Tsit5(), dtmax = 0.02, abstol = 1e-6, reltol = 1e-6)

# Modified version of Rozenweig-MacArthur model for predator-prey
# See https://arxiv.org/abs/2101.12107
function hassan_alkhayuon(u, p, t) # From talk of Hassan Alkhayuon
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
r1, r2 = 1.7, 2.5
# r, c, μ, ν, α, β, χ, δ = p
p = [r1, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]
ds = ContinuousDynamicalSystem(hassan_alkhayuon, u0, p)

# I will re-use code from InteractiveDynamics.jl, so I'll be defining some necessary
# quantities that may not make perfect sense without looking at the source code.
# Nevertheless the source code is quite clear and is just used to initialize some
# interactive observables and stuff.
tail = 500                  # length of plotted trajectory (in units of `dt`)
lims = ((0, 15), (0, 0.04)) # limits of state space
color = to_color("#1e1528") # trajectory color

figure = Figure(resolution = (1000, 800))
pinteg = DynamicalSystems.parallel_integrator(ds, [u0]; diffeq...)
obs, finalpoints = InteractiveDynamics.init_trajectory_observables(1, pinteg, tail, SVector(1,2), identity)

main = figure[1,1] = InteractiveDynamics.init_main_trajectory_plot(
    ds, figure, SVector(1,2), lims, pinteg, [color], obs, NamedTuple(), finalpoints, 1.0
)

# Plot the stability basin and attractor
active_color = to_color("#d7800a")
passive_color = RGBAf0(0.3, 0.4, 0.4, 0.5)
r1_color = Observable(active_color)
r2_color = Observable(passive_color)

function plot_limitcycle_basin!(main, r, color)
    set_parameter!(ds, 1, r)
    # plot unstable manifold of "e3" fixed point
    zz = SVector(p[3], 0.0)
    J = ds.jacobian(zz, p, 0)
    ev = eigvecs(J)
    # system in reverse time: just negative `f`.
    pprev = (u, p, t) -> -hassan_alkhayuon(u, p, t)
    dsrev = ContinuousDynamicalSystem(pprev, SVector(1.0, 1.0), p)
    e = ev[:, 1] # eigenvector correpsonding to positive eigenvalue
    tr = trajectory(dsrev, 3.5, zz .+ 1e-3 .* e; dtmax = 0.01)
    lines!(main, columns(tr)...; linewidth = 2.0, color = color, linestyle = :dash)
    # plot limit cycle
    tr = trajectory(ds, 20.0; dt = 0.01, Ttr = 100.0)
    lines!(main, columns(tr)...; linewidth = 4.0, color = color)
    return
end

plot_limitcycle_basin!(main, r2, r2_color)
plot_limitcycle_basin!(main, r1, r1_color)

# add controls for start/stop and toggle the r value below the figure
figure[2, 1] = controllayout = GridLayout(tellwidth = false)
runbutton = controllayout[1, 1] = Button(figure; label = "run")
rbutton = controllayout[1, 2] = Button(figure; label = "toggle r")
run, rtoggle = runbutton.clicks, rbutton.clicks

# Set up interactive part of the application: start/stop time evolution
isrunning = Observable(false)
on(run) do clicks; isrunning[] = !isrunning[]; end
on(run) do clicks
    @async while isrunning[]
        DynamicalSystems.step!(pinteg)
        ob = obs[1]
        last_state = DynamicalSystems.get_state(pinteg, 1)
        ob[] = push!(ob[], last_state) # push and trigger update with `=`
        finalpoints[] = [x[][end] for x in obs]
        sleep(0.00001) # use sleep to slow down animation
        isopen(figure.scene) || break # crucial, ensures computations stop if closed window
    end
end

# Set up interactive part of the application: toggle the r value
on(rtoggle) do clicks
    r = (r1, r2)[mod(clicks, 2)+1]
    set_parameter!(ds, 1, r)
    u_modified!(pinteg, true)
    x, y = r1_color[], r2_color[]
    r1_color[], r2_color[] = y, x
end

display(figure)
