# van_der_Pol_Poincare_section_v1.jl
# =============================================


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

# parameters and initial values
# -----------------------------
ntrans = 5000.   # transient periods
nattr = 2000    # no. of points on the attractor used for estimating winding no.
y_init = [0.1 0.1]

d = 5.
a = 1.
omega = 1.

p = [d a omega]
period = 2*pi / omega


 # transient to attractor
 # ----------------------
tspan = (0., ntrans*period)
prob = ODEProblem(van_der_Pol_ODE, y_init, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e12 )

# attractor
# ---------
tspan = (0., nattr*period)
y_init = sol.u[end]
prob = ODEProblem(van_der_Pol_ODE, y_init, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters = 1e12,  saveat = period)
  
u1vec = [u[1] for u in sol.u] 
u2vec = [u[2] for u in sol.u]   
nts = length(u1vec)

theta_vec = (atan.( u2vec , u1vec  ) .+ pi ) ./ (2*pi)  

# open plot
# ---------
# %%
fig = figure(figsize=(0.6*figx,figy)) 

ax1 = subplot(1,2,1)
plt.scatter(u1vec,u2vec,s=4, c = COLORS[1])

xticks([-2,0,2], ["-2","","2"])
yticks([-6,0,6], ["-6","","6"])
xlabel(L"u_1";labelpad = -20) 
ylabel(L"u_2";labelpad = -20) 

ax2 = subplot(1,2,2)
plt.scatter(theta_vec[1:end-1],theta_vec[2:end],s=4, c = COLORS[2])
xlabel(L"\theta_n"; labelpad = -20) 
ylabel(L"\theta_{n+1}"; labelpad = -20)
xticks([0,0.5,1], ["0", "", "2π"])
yticks([0,0.5,1], ["0", "", "2π"])
yticks([0,0.5,1])
xlim(0,1)
ylim(0,1)
add_identifiers!(fig)
fig.tight_layout(pad=0.33)

wsave(plotsdir("9", "vanderPol_poincare_section"), fig)
