# FitzHugh_Nagumo_2D_compute_spiral_wave_PS_v2.jl
# ================================================   u.p. 4.9.21

# Computes and saves spiral wave solution of the Fitzhugh-Nagumo model.
# Used for figures 11.4 and 11.6, but also in the animation spiralwave in animations/11.

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using Random
using Statistics

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

function FitzHugh_Nagumo_ODE(du,u,p,t)
    a, b, d, ε, hsq6 = p
    uu = @view u[1,:,:]
    vv = @view u[2,:,:]

    # Create padded u matrix to incorporate Newman boundary conditions - 9-point stencil
    uuu=[ [uu[2,2] uu[2,:]' uu[2,end-1] ] ; [ uu[:,2] uu uu[:,end-1] ] ; [ uu[end-1,2] uu[end-1,:]' uu[end-1,end-1] ] ]
    diff_term = d .* ( 
                4 .* (uuu[2:end-1,1:end-2] .+ uuu[2:end-1,3:end] .+ 
                        uuu[3:end,2:end-1] .+ uuu[1:end-2,2:end-1] ) .+ 
                        uuu[3:end,3:end] .+ uuu[3:end,1:end-2] .+ uuu[1:end-2,3:end] .+ 
                        uuu[1:end-2,1:end-2] .- 20 .*uu  )  ./ hsq6  

    du[1,:,:] = @. a*uu*(1 - uu)*(uu - b) - vv + diff_term
    du[2,:,:] = @. ε*(uu - vv)
end
    
L = 300.0  #  size of domain
N = 300   #  no. of grid points
h = L/N # spatial steps size
hsq6 = 6*h*h

# parameters
a = 3.    
b = 0.2   
d = 1.    #  diffusion
ε = 0.01
      
#  initial conditions
uu = 0.5 .* ones(N,N)
vv = 0.5 .* ones(N,N)
# This initial condition sets the spiral wave
uu[1:end, 1:(2N÷3)] .= 0
vv[1:(2N÷3 - N÷12), 1:end] .= 0

u0 = zeros(2,N,N)
u0[1,:,:] = uu
u0[2,:,:] = vv 

# solve ODE system
p = a, b, d, ε, hsq6 

saveat = 0.0:600.0
tspan = (0.0, saveat[end])

prob = ODEProblem(FitzHugh_Nagumo_ODE, u0, tspan, p)
@time sol = solve(prob, Tsit5();
    reltol=1e-6, abstol=1e-6, maxiters = Inf, saveat = saveat,
    progress = true, progress_steps = 10,
)
tvec = sol.t
uout = sol.u

params = @ntuple a b d ε hsq6 N L
data = @strdict tvec uout params 

wsave(datadir("FitzHugh", "spiralwave.jld2"), data)
