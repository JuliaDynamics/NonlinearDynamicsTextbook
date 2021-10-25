# Duffing_orbit_diagram_vs_period_single_v3.jl
# ============================================

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

function Poincare_orbit(d,a,omega,y_init,ttrans,nattr,pcol)

    p = [d a omega]
    period = 2*pi / omega
    ntrans = ceil( ttrans / period )   # no. of periods required for transient

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
    plt.scatter(period*ones(npts),u1vec,s=2, c = pcol)

    return y_init

end # of function



function plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)


color_line = COLORS[3]

omega  = 2. * pi / T_orbit
period = T_orbit

# transient to attractor
# ----------------------
ttrans = 100.
ntrans = ceil( ttrans / period )   # no. of periods required for transient
p = [d a omega]
tspan = (0., ntrans*period)
prob = ODEProblem(Duffing_ODE, y_init, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# orbit  
# -----
nsamp = 300
tspan = (0., nattr*period)
y_init = sol.u[end]
prob = ODEProblem(Duffing_ODE, y_init, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat = period/nsamp)

u1vec = [u[1] for u in sol.u]   
u2vec = [u[2] for u in sol.u]  
u1p = u1vec[end]

# plot line pointing to orbit
# ---------------------------
xx = [ T_orbit + 0.15delta_T  , T_orbit + linlen*delta_T ] 
yy = [ u1p + 0.15delta_u1p , u1p + linlen*delta_u1p ]
plt.plot(xx, yy, c = color_line, linewidth = 1)


# plot symbol indicating current solution
# ---------------------------------------
u1p = u1vec[end]
u2p = u2vec[end]
plt.scatter(T_orbit,u1p,s=55, c = color_symbol, zorder = 3 )

# plot orbit
# ----------
cc = [T_orbit + delta_T, u1p + delta_u1p] # center of orbit
# dd = [0.8 0.2]    # scaling

plt.plot(cc[1].+ dd[1] .* u1vec,  cc[2].+ dd[2] .* u2vec, c = color_orbit, linewidth=1.5)
end 


# ----------------------- main program -------------------------

@time begin

# parameters and initial values
# -----------------------------
d = 0.1
a = 7

T_min = 1 # 1
T_max = 29 # 32
n_T = 2000  # 1500
T_vec = [LinRange(T_min, T_max, n_T); LinRange(T_max, T_min, n_T)]

ttrans = 1000 # 150  # minimal transient time
nattr = 4     # no. of points on the attractor
y_init = [0.1 0.1]


# open plot
# ---------
fig = figure(figsize=(0.65*figx,figy))  

pcol = COLORS[1] 

#  control parameter loop (period)
#  ----------------------
for i_T = 1:length(T_vec)
    period = T_vec[i_T]
    omega = 2*pi / period
    y_init =  Poincare_orbit(d,a,omega,y_init,ttrans,nattr,pcol)
end # of period loop

# tight_layout()
ax = gca()
ax.set_ylim([-2.,1.7])
ax.set_xlim([T_min,T_max])


# Add orbits
# ==========
color_symbol = COLORS[4]
color_orbit = COLORS[4]

# Asymmetric orbit Ia
# ------------------
T_orbit = 6.
y_init = [0.5 -1.]
delta_T = 1.7
delta_u1p = 0.6 
linlen = 0.3
dd = [0.8 0.2]  * 0.6
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)


# Asymmetric orbit Ib
# ------------------
T_orbit = 6.
y_init = [-1. 1.]
delta_T = 1.2
delta_u1p = -0.8 
linlen = 0.6
dd = [0.8 0.2] * 0.6
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)


# Asymmetric orbit IIa
# --------------------
T_orbit = 15.
y_init = [1. 1.]
delta_T = -1.5
delta_u1p = 0.8 
linlen = 0.6
dd = [0.8 0.2]  
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)


# Asymmetric orbit IIb
# --------------------
T_orbit = 15.
y_init = [-1. 1.]
delta_T = -1.8
delta_u1p = -0.8 
linlen = 0.3
dd = [0.8 0.2]  
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)

# Symmetric orbit I
# -----------------
color_symbol = COLORS[6]
color_orbit = COLORS[6]

T_orbit = 17.
y_init = [-1. 1.]
delta_T = 1.
delta_u1p = 1. 
linlen = 0.5
dd = [0.8 0.2]  
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)


# Symmetric orbit II
# ------------------
T_orbit = 25.
y_init = [-1. 1.]
delta_T = 0.
delta_u1p = 0.8
linlen = 0.6
dd = [0.8 0.2] * 1.7
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)

# Symmetric orbit III
# ------------------
T_orbit = 3.
y_init = [-1. 1.]
delta_T = 0.
delta_u1p = 1.
linlen = 0.4
dd = [0.8 0.2] * 0.5
plot_orbit(a,d,T_orbit,y_init,delta_T,delta_u1p,color_symbol,color_orbit,linlen,dd)

end # of time 


xlabel(L"T = 2\pi / \omega")  
ylim(-2,2)
yticks([-2,2])
ylabel(L"u_1"; labelpad = -30)

fig.tight_layout(pad=0.3)

wsave(plotsdir("9", "duffing_orbitdiagram"), fig)
