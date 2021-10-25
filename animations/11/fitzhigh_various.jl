using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using Random, OrdinaryDiffEq
import GLMakie

function FitzHugh_Nagumo_ODE(du,u,p,t)
    a, b, d, epsilon, hsq6 = p
    uu = @view u[1,:,:]
    vv = @view u[2,:,:]

    # This code can be optimized a lot if time allows
    # Create padded u matrix to incorporate Newman boundary conditions - 9-point stencil
    uuu=[ [uu[2,2] uu[2,:]' uu[2,end-1] ] ; [ uu[:,2] uu uu[:,end-1] ] ; [ uu[end-1,2] uu[end-1,:]' uu[end-1,end-1] ] ]
    diff_term = d .* ( 
        4 .* (uuu[2:end-1,1:end-2] .+ uuu[2:end-1,3:end] .+ 
            uuu[3:end,2:end-1] .+ uuu[1:end-2,2:end-1] ) .+ 
            uuu[3:end,3:end] .+ uuu[3:end,1:end-2] .+ uuu[1:end-2,3:end] .+ 
            uuu[1:end-2,1:end-2] .- 20 .*uu  
    ) ./ hsq6  

    du[1,:,:] = a .* uu .*(1 .-uu).*(uu.-b) .- vv .+ diff_term
    du[2,:,:] = epsilon .* (uu .- vv )
end

# --------------------------  main program  -----------------------------

@time begin

# parameters
L = 300   # size of domain
N = 150   # no. of grid points
h = L/N # spatial steps size
hsq6 = 6*h*h

# FHN parameters
a = 3.    
b = 0.2   
d = 1.   #  diffusion
epsilon = 0.01 

# open figure
kplt = 0  # no. of subplot

allsols = []
allts = []

for ic = 1:2    # initial conditions ( = rows of the plot )

println("ic=",ic)

uu = zeros(N,N)
vv = zeros(N,N)   
Random.seed!(13371)
uu = uu .+ 0.01*randn(N,N)    # random perturbation

# local stimulius = concentric waves
if ic == 1 
   uu[40:41,40:41] .= 0.999
   uu[75:76,75:76] .= 0.999
end

# local stimulius = spiral waves
if ic == 2
   uu[1:end,10:11] .= 0.999
end   

# initial values   
u0 = zeros(2,N,N)
u0[1,:,:] = uu
u0[2,:,:] = vv 

p = a, b, d, epsilon, hsq6 

if ic == 1
    saveat = 0.0:2:600
    tspan = (saveat[1], saveat[end])

    prob = ODEProblem(FitzHugh_Nagumo_ODE, u0, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6, maxiters = 1e7, saveat)
    tvec = sol.t
    uout = [u[1, :, :] for u in sol.u]
    push!(allsols, uout)
    push!(allts, tvec)

elseif ic == 2
    saveat = 0.0:2.0:300.0
    tspan = (saveat[1], saveat[end])
    prob = ODEProblem(FitzHugh_Nagumo_ODE, u0, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6, maxiters = 1e7, saveat)

    uout = [u[1, :, :] for u in sol.u]
    # second perturbation line
    u = sol.u[end]
    u[1, 55:73,1:70] .= 0.999
    saveat = 300.0:2:900
    tspan = (saveat[1], saveat[end])
    prob = ODEProblem(FitzHugh_Nagumo_ODE, u, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6, maxiters = 1e7, saveat)

    tvec = sol.t
    uout2 = [u[1, :, :] for u in sol.u]
    append!(uout, uout2)
    push!(allsols, uout)
    push!(allts, tvec)

end
end # ic loop
end # cputime


# %% Animate
names = ("concentric", "plane_to_spiral")

for i in (1,2)
    name = names[i]
    uout = allsols[i]
    tvec = allts[i]

    heatobs = GLMakie.Observable(uout[1])
    tobs = GLMakie.Observable(0.0)
    titobs = GLMakie.lift(t -> "t = $(t)", tobs)

    fig = GLMakie.Figure(resolution = (600, 550))
    ax = fig[1,1] = GLMakie.Axis(fig; title = titobs)
    hmap = GLMakie.heatmap!(ax, heatobs; colormap = :tokyo, colorrange = (-0.2, 1))
    cb = GLMakie.Colorbar(fig[1, 2], hmap; width = 20)
    display(fig)

    GLMakie.record(
        fig, projectdir("animations", "11", "fitzhugh_$(name).mp4"), 
        1:length(tvec); framerate = 15
    ) do i
    
        tobs[] = tvec[i]
        heatobs[] = uout[i]
    end

end
