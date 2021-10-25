# Duffing resonance curves
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

function amplitude_curve(d,a,ωvec,ntrans,nattr,nres,y_init)
    nω = length(ωvec) # no. of driving frequencies (control parameter)
    amp_vec = zeros(nω)    # oscillation amplitudes (to be computed)
    ds = Systems.duffing(; d, f = a)

    for iω = 1:nω
        ω = ωvec[iω]
        set_parameter!(ds, 1, ω) # sets the ω parameter
        period = 2π/ω
        Δt = period/nres
        tr = trajectory(ds, nattr*period, y_init; Ttr = ntrans*period)
        amp_vec[iω] = maximum(u[1] for u in tr)
        y_init = tr[end]
    end
    return amp_vec
end

@time begin
# parameters and initial values
d = 0.1
omega_min = 0.02
omega_max = 3.3 
nomega = 1500 # 1500
ωvec = [LinRange(omega_min, omega_max, nomega); LinRange(omega_max, omega_min, nomega)]
avec = [1.1, 0.2, 0.04]
ccol = [COLORS[3], COLORS[2], COLORS[1]]
ntrans = 400   # no. of transient periods
nattr  = 2     # no. of periods on the attractor
nres   = 300   # resolution
amps = []

# compute resonance curves
for ia = 1:length(avec)
    a = avec[ia]
    y_init = [0.1, 0.1]
    amp_vec = amplitude_curve(d,a,ωvec,ntrans,nattr,nres,y_init)
    push!(amps, amp_vec)
end
end

# %% plot
fig = figure(figsize=(0.55*figx,figy))
ax = gca()
ax.set_xlabel(L"\omega") 
ax.set_ylabel("\$x_{\\mathrm{max}}\$")

for ia = 1:length(avec)
    a = avec[ia]
    ax.scatter(ωvec,amps[ia],s=2,c=ccol[ia], label = "\$a=$(avec[ia])\$")
end
ax.set_xlim([0,omega_max])
ax.legend(markerscale = 10)
fig.tight_layout(pad=0.3)

wsave(plotsdir("9", "duffing_resonance"), fig)
