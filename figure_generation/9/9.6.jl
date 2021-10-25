# van_der_Pol_orbit_diagram_vs_frequency_with_winding_number_v4.jl
# =================================================================

# with zoom
# new transient time

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using StaticArrays
using Random
using Statistics


function van_der_Pol_ODE(du,u,p,t)
    d, a, omega = p
    du[1] =  u[2]
    du[2] =  - u[1] - d*(u[1]*u[1]-1) * u[2] + a*sin(omega*t)
end


function Poincare_orbit(d,a,omega,y_init,ttrans,nattr,twind,pcol)

    p = [d a omega]
    period = 2*pi / omega
    ntrans = ceil( ttrans / period )   # no. of periods required for transient
    nwind = ceil( twind / period )   # no. of periods required for transient

 # transient to attractor
 # ----------------------
    tspan = (0., ntrans*period)
    prob = ODEProblem(van_der_Pol_ODE, y_init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e12 )

  # attractor  with root of Poincare map
  # ---------
    tspan = (0., nattr*period)
    y_init = sol.u[end]
    prob = ODEProblem(van_der_Pol_ODE, y_init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e12, saveat = period)

    u1vec = [u[1] for u in sol.u]   
    npts = length(u1vec)

    # Code for root of the Poincare map with 'saveat = period/2'
    # ----------------------------------------------------------
    # for j = 2:2:npts
    #     u1vec[j] = -u1vec[j]
    # end 

    # sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat = period/2)

    # u1vec = [u[1] for u in sol.u]   
    # npts = length(u1vec)
    # for j = 2:2:npts
    #     u1vec[j] = -u1vec[j]
    # end 


  # plot
  # ----
    plt.scatter(omega*ones(npts),u1vec,s=1, c = pcol)


  # compute winding number
  # ----------------------
    tspan = (0., nwind*period)
    # y_init = sol.u[end]
    prob = ODEProblem(van_der_Pol_ODE, y_init, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e12,  saveat = period/100)
  
    u1vec = [u[1] for u in sol.u] 
    u2vec = [u[2] for u in sol.u]   
    nts = length(u1vec)

    theta = 0.
    theta_old = atan.( u2vec[1] , u1vec[1])

    for its = 2:nts
        theta_new = atan.( u2vec[its] , u1vec[its]  )    # is dereasing = clockwise rotation
        # delta_th = theta_new  - theta_old  # < 0
        delta_theta = mod( theta_new  - theta_old  + pi, 2 .* pi) - pi
        theta = theta + delta_theta

        theta_old = theta_new
    end
    
    wind_no = - theta / ( 2 * pi * nwind )
    println(1 / wind_no)

    return y_init, wind_no

end # of function

# -------- main -----

println("this takes 56 min ..")

@time begin

# open plot
# ---------
fig = figure("van der Pol oscillator a=1",figsize=(figx,1.6*figy)) 

# parameters and initial values
# -----------------------------
d = 5.
a = 1.

ttrans = 6000   # minimal transient time
nattr = 500     # no. of points on the attractor
twind = 8000    # averaging of winding number
y_init = [0.1 0.1]


# compute and plot orbit diagram
# ------------------------------
omega_min = 0.2 
omega_max = 2. 
nomega =  1000 # 800 
omega_vec = LinRange(omega_min, omega_max, nomega)
wind_no_vec = zeros(nomega)  # winding numbers 

subplot(2,2,1)
   
pcol = COLORS[1] 
for iomega = 1:length(omega_vec)
    omega = omega_vec[iomega]
    y_init, wind_no =  Poincare_orbit(d,a,omega,y_init,ttrans,nattr,twind, pcol)
    wind_no_vec[iomega] = 1 / wind_no
end # of frequency loop

xlabel(L"\omega") 
ylabel(L"u_1") 

ax = gca()
ax.set_xlim([omega_min,omega_max])

# plot winding number
subplot(2,2,3)
plt.plot(omega_vec,wind_no_vec,c = "k")

xlabel(L"\omega") 
ylabel(L"W") 
ax = gca()
ax.set_xlim([omega_min,omega_max])

# zoom
# ====
omega_min = 0.8 
omega_max = 1. 
nomega =  1000 
omega_vec = LinRange(omega_min, omega_max, nomega)
wind_no_vec = zeros(nomega)  # winding numbers 

subplot(2,2,2)
   
pcol = COLORS[1] 
for iomega = 1:length(omega_vec)
    omega = omega_vec[iomega]
    y_init, wind_no =  Poincare_orbit(d,a,omega,y_init,ttrans,nattr,twind,pcol)
    wind_no_vec[iomega] = 1 / wind_no
end # of frequency loop

xlabel(L"\omega") 
ylabel(L"u_1") 

ax = gca()
ax.set_xlim([omega_min,omega_max])

# plot winding number
subplot(2,2,4)
plt.plot(omega_vec,wind_no_vec,c = "k")

xlabel(L"\omega") 
ylabel(L"W") 
ax = gca()
ax.set_xlim([omega_min,omega_max])

add_identifiers!(fig)
fig.tight_layout(pad=0.3)

end # of time 

# Output
# ------
# dname = "/Users/parlitz/ownCloud/NonlinearDynamicsTextbook/figure_generation/9/van_der_Pol/orbit_diagram/"
# fname = "van_der_Pol_orbit_diagram_vs_frequency_wind_no_v4_d=5_a=1_v5"  
# ftype = ".png"

# savefig(dname*fname*ftype)

