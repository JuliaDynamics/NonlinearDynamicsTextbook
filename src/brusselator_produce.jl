using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
using Random, OrdinaryDiffEq

# Run this file to integrate PDE Brusselator and save snapshots.
# Configure parameters HERE, including the `saveat` parameter
# (when to save the output fields). You can choose to call this file with a
# global existing variable `u0` (the 2xNxN array of the initial condition).
# if this is the case, change the `ic` parameter below to `:given`
a = 9.0
b = 10.2
L = 50.0 # size of domain
N = 250  # no. of grid points
d = 1.7
ic = :random
saveat = 1000:0.25:1025

# Remaining file (DO NOT ALTER!)
u0 = isdefined(Main, :u0) ? u0 : nothing
config = @strdict(a, b, d, L, N, ic, saveat, u0)

function Brusselator_ODE(du,u,p,t)
    a, b, d, hsq = p

    uu = @view u[1,:,:]
    vv = @view u[2,:,:]

    uE = circshift(uu, (0,1))   # u[:,[2:N 1]]
    uW = circshift(uu, (0,-1))  # u[:,[N 1:N-1]]
    uN = circshift(uu, (1,0))   # u[[N 1:N-1],:]
    uS = circshift(uu, (-1,0))  # u[[2:N 1],:]
    vE = circshift(vv, (0,1))
    vW = circshift(vv, (0,-1))
    vN = circshift(vv, (1,0))
    vS = circshift(vv, (-1,0))

    # finite difference formula
    au2v = @. a*uu*uu*vv
    du[1,:,:] = @. 1 - (b+1)*uu + au2v + (uE+uW+uN+uS - 4uu)/hsq 
    du[2,:,:] = @. b*uu - au2v + d*(vE+vW+vN+vS - 4vv)/hsq
    return
end

# This is for progress tracking of DiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

function produce_brusselator(config)
    @unpack a, b, d, saveat, N, L,ic = config
    h = L/N  # spatial step size
    hsq = h*h

    # initial condition
    if ic == :random
        us = 1  #  u steady state
        vs = b/a   # v steady state
        uu = fill(us, N, N)
        vv = fill(vs, N, N)
        Random.seed!(13371)
        uu = uu .+ 0.01*randn(N,N)    # random perturbation
        u0 = zeros(2,N,N)
        u0[1,:,:] .= uu
        u0[2,:,:] .= vv
    else
        u0 = config["u0"]
        @assert size(u0) == (2, N, N)
    end

    # Solve ODE
    p = (a, b, d, hsq)
    tspan = (0.0, saveat[end])
    @info "Initializing"
    prob = ODEProblem(Brusselator_ODE, u0, tspan, p)
    @info "Starting solve..."
    sol = solve(prob, Tsit5(); 
        reltol=1e-6, abstol=1e-6, maxiters = Inf, saveat,
        progress = true, progress_steps = 100,
    )
    t = sol.t
    u = sol.u
    output = @strdict u t a b d L N u0
    return output
end

prefix = length(saveat) > 10 ? "brusselator2D_dense" : "brusselator2D"

produce_or_load(datadir("Brusselator"), config, produce_brusselator;
    prefix, suffix = "jld2", tag = false, loadfile = false,
    force = true,
)
