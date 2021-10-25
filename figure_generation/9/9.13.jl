# Fig_GS_Roessler_Lorenz_v2.jl
# =============================

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using StaticArrays
using Random
using Statistics
using Neighborhood
using LinearAlgebra

function ccm(x, y, d, τ)
    ψ = copy(y); Mx = embed(x, d, τ); tw = Theiler(τ); 
    # bulk search of d+1 neighbors for all i points:
    idxs=bulkisearch(KDTree(Mx), Mx, NeighborNumber(d+1), tw)
    for i in 1:length(Mx)
      nni = idxs[i]; xi = Mx[i]; n1 = norm(xi - Mx[nni[1]])
      # implement equation (7.5):
      u = [exp(-norm(xi - Mx[j])/n1) for j ∈ nni]
      w = u ./ sum(u)
      ψ[i] = sum(w[k]*y[j] for (k, j) in enumerate(nni))
    end
    return cor(y, ψ) # Pearson correlation coef.
  end

function Roessler_Lorenz_aux_ODE(du,u,p,t)
    a, b, c, alpha,  cc = p

    du[1] = -u[2] - u[3] 
    du[2] = u[1] + a*u[2]      
    du[3] = b + u[3]*(u[1]-c)  

    
    du[1:3] = alpha .* du[1:3]   # adjust time scales
    
    du[4] = 10. * (-u[4] + u[5])  
    du[5] = 28*u[4] - u[5] - u[4]*u[6]  + cc*u[2]     
    du[6] = u[4]*u[5] - 2.666666 * u[6] 

    du[7] = 10. * (-u[7] + u[8]) 
    du[8] = 28*u[7] - u[8] - u[7]*u[9]  + cc*u[2]     
    du[9] = u[7]*u[8] - 2.666666 * u[9]     

end


function Roessler_Lorenz_CLE_ODE(du,u,p,t)
    a, b, c, alpha,  cc = p

    du[1] = -u[2] - u[3] 
    du[2] = u[1] + a*u[2]      
    du[3] = b + u[3]*(u[1]-c)  

     
    du[1:3] = alpha .* du[1:3]   # adjust time scales
    
    du[4] = 10. * (-u[4] + u[5]) 
    du[5] = 28*u[4] - u[5] - u[4]*u[6] + cc*u[2]
    du[6] = u[4]*u[5] - 2.666666 * u[6]

    # linearized equation
    du[7] =       -10 * u[7]  +   10 * u[8]
    du[8] = (28-u[6]) * u[7]         - u[8]     -  u[4] * u[9]    
    du[9] =      u[5] * u[7]  + u[4] * u[8]  - 2.666666 * u[9]     
    return nothing
end



#  ------------------ main program --------------------------


@time begin

# parameters and initial values
# -----------------------------
a = 0.2
b = 0.2  
c = 5.7  
alpha = 6 

ccmin = 0. 
ccmax = 12.
ncc = 100

ccvec = LinRange(ccmin,ccmax,ncc)
GSvec = zeros(ncc)
cxyvec = zeros(ncc)
cyxvec = zeros(ncc)

dxyvec = zeros(ncc)
dyxvec = zeros(ncc)

CLEvec = zeros(ncc)

ttrans = 2000  # transient time
tattr = 50000   # time for averaging on attractor

d = 6  
τ = 1  #  3
tsamp = 0.1 # 0.1 


# coupling parameter loop
# -----------------------
for icc = 1:ncc
    @show icc
    cc = ccvec[icc]
    p = a, b, c, alpha,  cc
    uinit = [-1.589656  3.0478997  0.823033  0.11  0.19  0.1 -0.1 -0.2 0.1]

   # transient to attractor
   # ----------------------
     tspan = (0., ttrans)
     prob = ODEProblem(Roessler_Lorenz_aux_ODE, uinit, tspan, p)
     sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e7)

   # attractor  
   # ---------
 
    tspan = (0., tattr)
    y_init = sol.u[end]
    prob = ODEProblem(Roessler_Lorenz_aux_ODE, y_init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e7, saveat = tsamp)

    tvec = sol.t
    ures  = sol.u
    nts = length(tvec)

    xpts = zeros(3,nts)
    ypts = zeros(3,nts)
    zpts = zeros(3,nts)

    for its = 1:nts     # .... this can be done more effciently
        uu = ures[its]
        xpts[1:3,its] = uu[1:3]
        ypts[1:3,its] = uu[4:6] 
        zpts[1:3,its] = uu[7:9]
    end

   # Synchronization error (auxiliary system method)
   # ----------------------------------------------- 
    GSerror = 0.
    for its = 1:nts # .... can this be done without a loop?
        GSerror = GSerror +  norm(ypts[:,its] - zpts[:,its] )
    end
    GSerror = GSerror / nts
    GSvec[icc] = GSerror

    # convergent cross mapping
    # ------------------------
    xx = xpts[1,:]   
    yy = ypts[1,:]

    cxyvec[icc] = ccm(xx, yy, d, τ)
    cyxvec[icc] = ccm(yy, xx, d, τ)

    if icc == 1
    fig = plt.figure("Fig Roessler Lorenz d=3 τ=1",figsize=(figx,figy))

    subplot(1,2,1)
    plt.scatter(xx[1:end-1],xx[2:end],s=1)
    
    subplot(1,2,2)
    plt.scatter(yy[1:end-1],yy[2:end],s=1)
    end


   # CLEs
   # ----
   tsamp = 1. # 0.1 
   tspan = (0., tattr)
   y_init = sol.u[end]
   y_init[7:9] = [1 0 0]
   ncle = 50000 # Lyapunov steps
   csum = 0.

   for icle = 1:ncle
       tspan = (0., tsamp)
       prob = ODEProblem(Roessler_Lorenz_CLE_ODE, y_init, tspan, p)
       sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e7, saveat = tsamp)

       tvec = sol.t
       y_init = sol.u[end]

       ynorm = sqrt( norm( y_init[7:9] ) ) 
       y_init[7:9] = y_init[7:9] ./ ynorm   # normalization of tangent vector

       if icle > 300
          csum = csum + log(ynorm)
       end

   nts = length(tvec)

   end

   CLE = csum / (ncle-300)
   CLEvec[icc] = CLE

# println("c=",cc,"  GS=",GSerror, "  ccm_xy=", cxyvec[icc],"  ccm_yx=",cyxvec[icc],"  dxy=",dxyvec[icc],"  dyx=",dyxvec[icc],"  CLE=", CLE)
println("c=",cc,"  GS=",GSerror, "  ccm_xy=", cxyvec[icc],"  ccm_yx=",cyxvec[icc],"  CLE=", CLE)



end


end  # of time

# %% plot

fig = plt.figure(figsize=(0.5*figx,figy))

plt.plot(ccvec,GSvec ./maximum(GSvec),color=COLORS[1])  # ,linewidth=1)
plt.plot(ccvec,CLEvec ./ maximum(CLEvec),color=COLORS[2], ls = "--") 

plt.plot(ccvec,cxyvec,color=COLORS[3])
plt.plot(ccvec,cyxvec,color=COLORS[4], ls = "-.")

plt.xlabel(L"k"; labelpad = -20) 


ax = gca()
ax.set_xlim([ccmin,ccmax])
ax.set_xticks(1:2:11)
# ax.set_ylim([-2.1,2.6])
ax.legend([L"E"; L"\lambda_1^{\rm{CLE}}";
           "CCM: \$x_1 \\to y_1\$"; "CCM: \$y_1 \\to x_1\$"],
         bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
        ncol=2, mode="expand", borderaxespad=0, handlelength=2,
)

fig.tight_layout(pad=0.3)
wsave(plotsdir("9", "generalized_sync"), fig)
