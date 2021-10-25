# vanderPol orbit diagram and winding number versus ω
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

function stroboscopic_plot_winding(ds,ω,y_init,Ttr,nattr,nwind)
    period = 2π/ω
    ds.p[3] = period
    T = period*nattr   # no. of periods required for transient
    Ttr = period*ceil(Ttr/period)     # no. of periods used for averaging

    tr = trajectory(ds, T, y_init; Δt = period, Ttr)
    u1ret = tr[:, 1]

    # compute winding number
    Δt = period/20
    tr = trajectory(ds, nwind*period, y_init; Δt, Ttr)
    u1, u2 = columns(tr)
    θ = 0.0
    θ_old = atan(u2[1], u1[1])
    for (x, y) ∈ zip(u1, u2)
        θ_new = atan(y, x)
        θ += mod(θ_new - θ_old + π, 2π) - π
        θ_old = θ_new
    end
    W = abs(θ / nwind*2π)
    return tr[end], W, u1ret
end # of function

# -------- main -----

d = 5.0
a = 1.0
Ttr = 1000      # minimal transient time
nattr = 1000    # no. of points on the attractor
nwind = 100    # averaging of winding number
y_init = SVector(0.1, 0.1)
ds = Systems.vanderpol(y_init; μ = d, F = a)


# compute and plot orbit diagram
# ------------------------------
ω_min = 0.2 
ω_max = 2. 
nω =  100 # 800 
ω_vec = LinRange(ω_min, ω_max, nω)
wind_no_vec = zeros(nω)  # winding numbers 
u1s = []
   
pcol = COLORS[1] 
for iω = 1:length(ω_vec)
    @show iω
    ω = ω_vec[iω]
    y_init, wind_no, u1 =  stroboscopic_plot_winding(ds,ω,y_init,Ttr,nattr,nwind)
    wind_no_vec[iω] = 1 / wind_no
    push!(u1s, u1)
end # of frequency loop

# zoom
ω_min = 0.8 
ω_max = 1.0
ω_vec_zoom = LinRange(ω_min, ω_max, nω)
wind_no_vec_zoom = zeros(nω)  # winding numbers 
u1s_zoom = []


for iω = 1:length(ω_vec_zoom)
    @show iω
    ω = ω_vec_zoom[iω]
    y_init, wind_no, u1 =  stroboscopic_plot_winding(ds,ω,y_init,Ttr,nattr,nwind)
    wind_no_vec_zoom[iω] = 1 / wind_no
    push!(u1s_zoom, u1)
end # of frequency loop


# %% plot
fig, axs = subplots(2,2; figsize = (figx, 1.5figy))

j = 0
for (uvec, ωvec) in zip((u1s, u1s_zoom), (ω_vec, ω_vec_zoom))
    j += 1
    for (i, ω) in enumerate(ωvec)
        axs[1,j].scatter(fill(ω, length(uvec[i])), uvec[i], s = 1, c = "C0", alpha = 0.2)
    end
    w = (wind_no_vec, wind_no_vec_zoom)[j]
    axs[2,j].plot(ω_vec, w; color = "C2")
    axs[1,j].set_xticklabels([])
    axs[1,j].set_xlim(ωvec[1], ωvec[end])
    axs[2,j].set_xlim(ωvec[1], ωvec[end])
end






axs[1,1].

subplot(2,2,1)
xlabel(L"\ω") 
ylabel(L"u_1") 

ax = gca()
ax.set_xlim([ω_min,ω_max])

# plot winding number
subplot(2,2,3)
plt.plot(ω_vec,wind_no_vec,c = "k")

xlabel(L"\ω") 
ylabel(L"W") 
ax = gca()
ax.set_xlim([ω_min,ω_max])


xlabel(L"\ω") 
ylabel(L"u_1") 

ax = gca()
ax.set_xlim([ω_min,ω_max])

# plot winding number
subplot(2,2,4)
plt.plot(ω_vec,wind_no_vec,c = "k")

xlabel(L"\ω") 
ylabel(L"W") 
ax = gca()
ax.set_xlim([ω_min,ω_max])

add_identifiers!(fig)
fig.tight_layout(pad=0.3)


# Output
# ------
# dname = "/Users/parlitz/ownCloud/NonlinearDynamicsTextbook/figure_generation/9/van_der_Pol/orbit_diagram/"
# fname = "van_der_Pol_orbit_diagram_vs_frequency_wind_no_v4_d=5_a=1_v5"  
# ftype = ".png"

# savefig(dname*fname*ftype)

