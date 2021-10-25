using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random, OrdinaryDiffEq

a = 0.2
b = 0.2
c = 5.7
amu = 0.02
ω1 = 1. + amu
ω2 = 1. - amu
kl = 201
kmin = 0
kmax = 0.05

config = @strdict ω1 ω2 a b c kmin kl kmax

function produce_croesslers_spectrum(config)

@unpack ω1, ω2, a, b, c, kmin, kl, kmax = config

y_init = [0.0, 0.2, 0., 0.11, 0.19, 0.1]

croesslers = Systems.coupled_roessler(y_init; ω1, ω2, a, b, c)

diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9, maxiters = Inf)

Δt = 0.01
T = 1000
t = 0:Δt:T

krange = range(kmin, kmax; length = kl)

λs = zeros(length(krange), 4)
sync_error= zeros(length(krange))

for (i, k) in enumerate(krange)
    @show i, k
    set_parameter!(croesslers, 6, k)
    set_parameter!(croesslers, 7, k)
    λs[i, :] .= lyapunovspectrum(croesslers, 100000, 4; Ttr = 5000, diffeq...)

    # computation of phase difference -----
    T = 20000 # length of time series 
    Δt = 0.1
    X = trajectory(croesslers, T; Ttr = 5000, Δt, diffeq...)  # which initial conditions are used?
    t = 0:Δt:T
    x1vec = X[:, 1]
    x2vec = X[:, 2]
    x4vec = X[:, 4]
    x5vec = X[:, 5]

    theta1 = 0.
    theta1_old = atan.( x2vec[1] , x1vec[1] )
    theta2 = 0.
    theta2_old = atan.( x5vec[1] , x4vec[1] )
    
    nts = length(t)
    for its = 2:nts
        theta1_new = atan.( x2vec[its] , x1vec[its]  )    # is decreasing = clockwise rotation
        delta_theta1 = mod( theta1_new  - theta1_old  + pi, 2 .* pi) - pi

        theta1 = theta1 + delta_theta1
        theta1_old = theta1_new
 
        theta2_new = atan.( x5vec[its] , x4vec[its]  )    # is decreasing = clockwise rotation
        delta_theta2 = mod( theta2_new  - theta2_old  + pi, 2 .* pi) - pi
        theta2 = theta2 + delta_theta2
        theta2_old = theta2_new
    end
    
    sync_error[i]  = abs(theta1 - theta2) / (t[end] - t[1])
end

ret = @strdict λs sync_error
ret = merge(config, ret)
return ret
end

file, path = produce_or_load(
    datadir("croesslers"), config, produce_croesslers_spectrum; 
    prefix = "croesslers_spectrum", tag = false,
)

@unpack λs, sync_error = file


# %% Plot
fig = figure(figsize = (0.66figx, 1.2figy))
gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[2,1])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])

ax1.plot(krange, λs)
# ax.set_xlabel("\$k\$")        
ax1.set_xticklabels([])
ax1.set_ylabel("\$\\lambda_i\$")
ax1.set_xlim(kmin,kmax)
ax1.set_ylim(-0.07,0.08)

ax2.plot(krange, sync_error)
ax2.set_xlabel("\$k\$"; labelpad = -20)
ax2.set_ylabel("\$\\Delta \\Omega\$")
ax2.set_xlim(krange[1],krange[end])

add_identifiers!(fig)

fig.tight_layout(;pad = 0.33)
wsave(plotsdir("9", "coupled_roessler_spectrum"), fig)
