# Chimera_states_v4.jl
# ====================   u.p. 13.9.21

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

function chimera_ODE(du,u,p,t)
    ω, R, alpha, nosc = p
    @inbounds for i = 1:nosc
        term = 0.0
        for j = i-R:i+R
            jj = 1 + mod(j-1,nosc)
            term = term + sin( u[i] - u[jj] + alpha )
        end
        du[i] = ω - term / (2*R )
    end
    return
end


#  ------------------------------  main program  --------------------------------
fig, axs = subplots(1,3; figsize=(figx,figy))

for axi in 1:2
if axi == 2
    # case a:  with chaotic transient
    nosc = 40 # 50#  40 # 100 # 40    # number of oscillators 
    R = 14 # 15 # 14 # 35  # 14  # coupling range
    ish = 102  # select snapshot
    tsamp = 50.  # sampling time
    tspan = 25000.  # total lengthof time series
else
    # case b:  no chaotic transient (or transient end much later), but more oscillators
    nosc = 100 # 50#  40 # 100 # 40    # number of oscillators 
    R = 35 # 15 # 14 # 35  # 14  # coupling range
    ish = 102  # select snapshot
    tsamp = 50.  # sampling time
    tspan = 25000.  # total lengthof time series
end

indx = LinRange(1,nosc,nosc)
ω = 0.
alpha = 1.46

p = ω, R, alpha, nosc

iseed = 7117 # 5117 
rng = MersenneTwister(iseed) # random initial conditiosn in [0 ,2*pi]
u0 = 2 .* pi .* (rand(rng,nosc) .- 0.5)

prob = ODEProblem(chimera_ODE, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat = tsamp)
uout = sol.u
tvec = sol.t
nts = length(tvec)

if axi == 1

axs[1].scatter(indx,mod.(uout[ish], 2*pi),color=COLORS[1])   
axs[1].set_xlim(indx[1],indx[end])
axs[1].set_xticks(indx[1]:30:indx[end])
axs[1].set_ylim(0,2pi)
axs[1].set_xlabel(L"i"; labelpad = -20) 
axs[1].set_ylabel(L"\phi"; labelpad=-10) 
end

# plot evolution of phase velocity
# compute phase velocities (could perhaps be replaced by call of ODE function)
qmat = zeros(nosc,nts) 
for its = 1:nts
    u = mod.(uout[its], 2*pi)
    for i = 1:nosc
        term = 0.0
        for j = i-R:i+R
            jj = 1 + mod(j-1,nosc)
            term = term + sin( u[i] - u[jj] + alpha )
        end
        qmat[i,its] = ω - term / (2*R )
    end
end

im1 = axs[1+axi].pcolormesh(tvec,indx,qmat,cmap = "inferno") 
axs[1+axi].set_xlabel(L"t"; labelpad = -20) 
axs[axi+1].set_ylabel(L"i"; labelpad = -20)
axs[axi+1].set_yticks([indx[1], indx[end]])
axs[axi+1].set_xticks([0, tvec[end]])


if axi == 2
    cbar = plt.colorbar(im1)  # , ticks=[-0.2, 0.2, 0.6, 1.], orientation="horizontal")
    # cbar.set_clim(1., 4.) 
    cbar.set_label(L"d \phi /  dt") 
end

end

# %% 
add_identifiers!(fig, axs; xloc = 0.05)

fig.tight_layout(pad=0.3)

wsave(plotsdir("10", "chaimera"), fig)
