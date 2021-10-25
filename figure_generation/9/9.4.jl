# Autonomous van der Pol
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

# parameters and initial values
dvec = [0.1, 5]    # damping parameters
y_0 = [0.1, 0.1]   # initial condition
ttrans = 100   # transient time
tattr = 16     # time on attractor
Δt = 0.01   # sampling time 
tvec = 0:Δt:tattr

fig, axs = subplots(1,2; figsize=(0.6figx,0.8figy)) 

for id = 1:2
    d = dvec[id]
    y_init = y_0
    ds = Systems.vanderpol(; F = 0, μ = d)
    tr = trajectory(ds, tattr, y_0; Δt, Ttr = ttrans)
    x1vec, x2vec = columns(tr)
    axs[1].plot(tvec,x1vec; c = COLORS[id] ) 
    axs[2].plot(x1vec,x2vec; c = COLORS[id] )  
end
axs[1].set_xlabel(L"$t$"; labelpad = -20) 
axs[1].set_ylabel(L"$x$", rotation = 0) 
axs[2].set_xlabel(L"$x$"; labelpad = -20)  
axs[2].set_ylabel(L"$\dot x$", rotation = 0)
axs[1].set_xticks([0,15])
axs[2].set_xticks([-2,0,2])
axs[2].set_xticklabels(["-2", "", "2"])
add_identifiers!(fig)
fig.tight_layout(pad=0.33)

wsave(plotsdir("9", "vanderPol_autonomous"), fig)
