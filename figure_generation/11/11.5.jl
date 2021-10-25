# FitzHugh_Nagumo_2D_concentric_plane_spiral_waves_v1.jl
# ======================================================   u.p. 4.8.21

# only concetric and spiral waves - no linear wave

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using Random

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
   tspan = (0., 300.)
   tsamp = [0., 90., 150., 300.]
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
   prob = ODEProblem(FitzHugh_Nagumo_ODE, u0, tspan, p)
   sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, maxiters = 1e7, saveat = tsamp)
   tvec = sol.t
   uout = sol.u
   push!(allsols, uout)
   push!(allts, tvec)

elseif ic == 2
   tspan = (0.,300.)
   tsamp = [300.]
 
   prob = ODEProblem(FitzHugh_Nagumo_ODE, u0, tspan, p)
   sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, maxiters = 1e7, saveat = tsamp)
 
   # second perturbation line
   uout = sol.u
   u = uout[end]
   u[1,55:73,1:70] .= 0.999
   tspan = (300.,1000.)
   tsamp = [300., 400., 780., 781.]
   prob = ODEProblem(FitzHugh_Nagumo_ODE, u, tspan, p)
   sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, maxiters = 1e7, saveat = tsamp)

   tvec = sol.t
   uout = sol.u

   tvec[2:4] = tvec[1:3]
   tvec[1] = 0.

   uout[2:4] = uout[1:3]
   uout[1] = u0  # initial condition
   push!(allsols, uout)
   push!(allts, tvec)

end
end # ic loop
end # cputime

fig = plt.figure(figsize=(figx,1.58*figy))

kplt = 0
for ic in (1,2)
   uout = allsols[ic]
   tvec = allts[ic]
   nts = length(tvec)

   for n = 1:nts
      ures = uout[n]
      kplt = kplt + 1
      subplot(2,4, kplt)
      im = plt.pcolormesh(ures[1,:,:],cmap = "turbo")   # "gnuplot" ) #"inferno")
      im.set_clim(-0.2, 1.)

      ax = gca()
      ax.set_xticks([]) 
      ax.set_yticks([]) 
      plt.title(L"t = "*string(floor(Int,tvec[n])), fontsize = 24)
      ax.set_aspect("equal")
   end

end # end of ic-loop

fig.tight_layout(pad=0.3)
add_identifiers!(fig)

wsave(plotsdir("11", "FitzHugh_Nagumo_2D_concentric_plane_spiral_waves"), fig)
