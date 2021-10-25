using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
using Random, OrdinaryDiffEq

# Run this file to integrate PDE Brusselator and save snapshots.
# Configure parameters HERE, including the `saveat` parameter
# (when to save the output fields). You can choose to call this file with a
# global existing variable `u0` (the 2xNxN array of the initial condition).
# if this is the case, change the `ic` parameter below to `:given`
a = 3.0
b = 0.2
ε = 0.01
L = 300.0 # size of domain
N = 600  # no. of grid points
d = 1.0
u0 = :spiralwave
saveat = 1000:0.25:1025

# Remaining file (DO NOT ALTER!)
config = @strdict(a, b, d, L, N, saveat, u0)

function FitzHugh_Nagumo_ODE(du,u,p,t)
    a, b, d, hsq6 = p

    uu = @view u[1,:,:]
    vv = @view u[2,:,:]

    # TODO: What is this crazyiness? Creating `uuu` costs as much as computing
    # the trajectories of the solar system. But why though?
    # Aren't the Newman boundary conditions literally setting du[1, 1, 1] = du[1, N, N] = 0?

    # Create padded u matrix to incorporate Newman boundary conditions - 9-point stencil
    uuu=[ [uu[2,2] uu[2,:]' uu[2,end-1] ] ; [ uu[:,2] uu uu[:,end-1] ] ; [ uu[end-1,2] uu[end-1,:]' uu[end-1,end-1] ] ]
    diff_term = d .* ( 
                4 .* (uuu[2:end-1,1:end-2] .+ uuu[2:end-1,3:end] .+ 
                        uuu[3:end,2:end-1] .+ uuu[1:end-2,2:end-1] ) .+ 
                        uuu[3:end,3:end] .+ uuu[3:end,1:end-2] .+ uuu[1:end-2,3:end] .+ 
                        uuu[1:end-2,1:end-2] .- 20 .*uu  )  ./ hsq6  


    # finite difference formula
    @. du[1,:,:] = a*uu*(1 - uu)*(uu - b) - vv + diff_term    
    @. du[2,:,:] = ε*(uu - vv)
    return
end

# This is for progress tracking of DiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

function produce_brusselator(config)
    @unpack a, b, ε, d, saveat, N, L, u0 = config
    h = L/N  # spatial step size
    hsq6 = 6*h*h

    # initial condition
    if u0 == :spiralwave
        uu = fill(0.5, N, N)
        vv = fill(0.5, N, N)
        uu[1:end, 1:(2N÷3)] .= 0
        vv[1:(2N÷3 - N÷12), 1:end] .= 0
        u0 = zeros(2,N,N)
        u0[1,:,:] = uu
        u0[2,:,:] = vv 
    end
    @assert u0 isa AbstractArray
    # Solve ODE
    p = a, b, d, epsilon, hsq6 
    tspan = (0.0, maximum(saveat))
    @info "Initializing"
    prob = ODEProblem(FitzHugh_Nagumo_ODE, u0, tspan, p)
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

prefix = length(saveat) > 10 ? "fitzhugh_dense" : "fitzhugh"

produce_or_load(datadir("FitzHugh"), config, produce_brusselator;
    prefix, suffix = "jld2", tag = false, loadfile = false,
    force = true,
)
