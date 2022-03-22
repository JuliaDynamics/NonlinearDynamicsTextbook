# SIR_v2.jl
# ===============================
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
#using DynamicalSystems 
using PyPlot, OrdinaryDiffEq
using StaticArrays


function SIR_ODE(u,p,t)
    beta, gamma, N = p
    S, I, R = u
    dS =  - beta * I * S / N 
    dI =  beta * I * S / N  - gamma * I
    dR = gamma * I
    return SVector(dS, dI, dR)
end

N = 10000
beta = 0.3
gamma = 0.1 
p = beta, gamma, N

u0 = SVector(N-1., 1., 0.)

tspan = (0., 120)
    
prob = ODEProblem(SIR_ODE, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,tsamp = 0.1)

tvec = sol.t
Svec = [u[1] for u in sol.u]
Ivec = [u[2] for u in sol.u]
Rvec = [u[3] for u in sol.u]

# plot
# ----
fig = figure("SIR",figsize=(0.5*figx,figy)) 

plot(tvec,Svec/N, c = COLORS[1] ) 
plot(tvec,Rvec/N, c = COLORS[3] ) 
plot(tvec,Ivec/N, c = COLORS[2] ) 

xlabel(L"$t$") 
ax = gca()
ax.set_xlim([0.,120.])
# ax.set_ylim([-2.1,2.5])
ax.legend([L"S/N"; L"R/N"; L"I/N";], fontsize = 25)
fig.tight_layout(pad=0.3)

# save
# ----
# dname = "/Users/parlitz/ownCloud/NonlinearDynamicsTextbook/figure_generation/10/SIR/"
# fname = "SIR_model"
# ftype = ".png"

wsave(plotsdir("10", "SIR_model"), fig)