# Aliev_Panfilov_2D_compute_v4.jl
# ================================   u.p. 5.8.21  / 5.9.21

# Aliev-Panfilov model
# with ODE-solver 
# compute a save solution

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using PyPlot, OrdinaryDiffEq
using JLD
using Random
using Statistics
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

function Aliev_Panfilov_ODE(du,u,p,t)
    d, k, a, epsilon, mu1, mu2, hsq6  = p

    uu = u[1,:,:]
    vv = u[2,:,:]
    uuu=[ [uu[2,2] uu[2,:]' uu[2,end-1] ] ; [ uu[:,2] uu uu[:,end-1] ] ; [ uu[end-1,2] uu[end-1,:]' uu[end-1,end-1] ] ]

    # 9-point stencil
    diff_term = d .* ( 
                4 .* (uuu[2:end-1,1:end-2] .+ uuu[2:end-1,3:end] .+ 
                      uuu[3:end,2:end-1] .+ uuu[1:end-2,2:end-1] ) .+ 
                    uuu[3:end,3:end] .+ uuu[3:end,1:end-2] .+ uuu[1:end-2,3:end] .+ 
                    uuu[1:end-2,1:end-2] .- 20 .*uu  )  ./ hsq6  
    
    du[1,:,:] = k .* uu .*(1 .-uu).*(uu.-a) .- uu.*vv .+ diff_term 
    du[2,:,:] = (epsilon .+ (mu1 .* vv) ./ (mu2 .+ uu))  .* (.- vv .- k .* uu .* (uu .- a .- 1 ) )
end

@time begin


# new simulation starting from initial conditions
# -----------------------------------------------
L = 100  #  size of domain
N = 400  #  no. of grid points
h = L/(N-1) # spatial steps size
hsq6 = 6*h*h

# parameters
k = 8
a = 0.05    
mu1 = 0.2   
mu2 = 0.3
d = 0.2   #  diffusion
epsilon = 0.002 

# initial values
# ---------------
uu = zeros(N,N)  #  u steady state
vv = zeros(N,N)  #  v steady state
uu[1:end,1:20] .= 0.9   # N = 400
# uu[1:end,1:10] .= 0.9   # N = 200


u0 = zeros(2,N,N)
u0[1,:,:] = uu
u0[2,:,:] = vv 

p = d, k, a, epsilon, mu1, mu2, hsq6 
saveat = (0.0, 120.0) # [0. , 20. , 120.]
tspan = (0., saveat[end])

prob = ODEProblem(Aliev_Panfilov_ODE, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat = saveat,
progress = true, progress_steps = 100,
)
uout = sol.u

# perturbation
# -------------
ures = uout[end]
uu = ures[1,:,:]
vv = ures[2,:,:]

# uu[1:30,:] .= 1.  # N = 200
uu[1:60,:] .= 1.  # N = 400

u0[1,:,:] = uu
u0[2,:,:] = vv   

# evolve forward
# --------------
# saveat = [700. , 710., 720., 730., 740., 750., 760., 770., 780.]
# tspan = (0., saveat[end])
saveat = [600]
tspan = (0., saveat[end])
    
prob = ODEProblem(Aliev_Panfilov_ODE, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat = saveat,
progress = true, progress_steps = 100,
)

# compute snapshots
# -----------------
uout = sol.u
ures = uout[end]

saveat = 1.0
tspan = (0., 200)
    
prob = ODEProblem(Aliev_Panfilov_ODE, ures, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat = saveat,
progress = true, progress_steps = 100,
)


end # @time

tvec = sol.t
uout = sol.u
nts = length(tvec)

# plot
# ----
fig = plt.figure("AP_2D",figsize=(figx,2*figy))
iplt = 0

for n = 1:9
    iplt = iplt + 1
    ures = uout[n]
    subplot(3,3,iplt)
    im = plt.pcolormesh(ures[1,:,:],cmap = "inferno")  # "gnuplot" ) #"inferno")

    ax = gca()
    ax.set_xticks([]) 
    ax.set_yticks([]) 

    plt.title(L"t = "*string(floor(Int,tvec[n])), fontsize = 24)
end

fig.tight_layout(pad=0.3)
path = datadir("Aliev_Panfilov_snapshots.jld")
# save(path, "d", d, "k", k, "a", a, "epsilon", epsilon, "mu1", mu1, "mu2", mu2, "L", L, "N", N, "tvec", tvec, "uout", uout)

