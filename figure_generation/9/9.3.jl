# Duffing_orbit_diagram_vs_period_zooom_v1.jl
# ===========================================

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using StaticArrays
using Random
using Statistics

function Duffing_ODE(du,u,p,t)
    d, a, omega = p
    du[1] =  u[2]
    du[2] =  - u[1] - u[1]*u[1]*u[1] - d*u[2] + a*sin(omega*t)
end

function Poincare_orbit(d,a,omega,y_init,ntrans,nattr,pcol)

    p = [d a omega]
    period = 2*pi / omega

 # transient to attractor
 # ----------------------
    tspan = (0., ntrans*period)
    prob = ODEProblem(Duffing_ODE, y_init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

  # attractor  with root of Poincare map
  # ---------
    tspan = (0., nattr*period)
    y_init = sol.u[end]
    prob = ODEProblem(Duffing_ODE, y_init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat = period/2)

    u1vec = [u[1] for u in sol.u]   
    npts = length(u1vec)
    for j = 2:2:npts
        u1vec[j] = -u1vec[j]
    end 

  # plot
  # ----
    plt.scatter(period*ones(npts),u1vec,s=1, c = pcol)

    return y_init

end # of function


# parameters and initial values
# -----------------------------
d = 0.1
avec = [19.7 19.75]
n_T = 1500

ntrans = 800   # no. of transient periods
nattr = 64     # no. of points on the attractor
y_init = [1 0]


# open plot
# ---------
fig = figure() 

# subplot loop
# ------------
@time for iplt = 1:2

    if iplt  == 1
        a = avec[1]
        T_min = 4.1 # 1
        T_max = 5.8 # 32
        T_vec = LinRange(T_min, T_max, n_T) 
    else
        a = avec[2]
        T_min = 5.06 # 1
        T_max = 5.42 # 32
        T_vec = LinRange(T_max, T_min, n_T) 
    end

    subplot(1,2,iplt)

    xlabel(L"T = 2\pi / \omega") 
    ylabel(L"u_1") 
    pcol = COLORS[1] 

#  control parameter loop (period)
#  ----------------------
   for i_T = 1:length(T_vec)
       period = T_vec[i_T]
       omega = 2*pi / period
       y_init =  Poincare_orbit(d,a,omega,y_init,ntrans,nattr,pcol)
   end # of period loop

#   tight_layout()
#   ax = gca()
Tbmin = 5.1
Tbmax = 5.3
u1bmin = -0.94
u1bmax = -0.7 


#  plot zoom box
if iplt == 1
#    xbox = [Tbmin, Tbmax, Tbmax , Tbmin, Tbmin]
#    ybox = [u1bmin, u1bmin, u1bmax, u1bmax, u1bmin]
#    plt.plot(xbox,ybox,linewidth=1.0,color="k")
   ax = gca()
   ax.set_xlim([T_min,T_max])
   ax.set_ylim([-1.7,1.4])
end
if iplt == 2
    ax = gca()
    ax.set_xlim([Tbmin,Tbmax])
    ax.set_ylim([u1bmin,u1bmax])
end


end # of subplot loop

# %%

add_identifiers!(fig)
fig.tight_layout(pad=0.3)

wsave(plotsdir("9", "duffing_orbitdiagram_zoom"), fig)
