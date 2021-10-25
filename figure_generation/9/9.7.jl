# vanderPol winding number arnold tongues
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

function winding_number(ds, ω, a, u0, Ttr, T)
    period = 2π/ω
    ds.p[2:3] .= (a, period)
    # T = period*ceil(T/period)   # no. of periods required for transient
    # Ttr = period*ceil(Ttr/period)     # no. of periods used for averaging
    tr = trajectory(ds, T, u0; Δt = period/20, Ttr)
    tvec = 0:period/20:T
    u1, u2 = columns(tr)
    θ = 0.0
    θ_old = atan(u2[1], u1[1])
    for (x, y) ∈ zip(u1, u2)
        θ_new = atan(y, x)
        θ += mod(θ_new - θ_old + π, 2π) - π
        θ_old = θ_new
    end
    W = abs(θ / (tvec[end] * 2π/period))
    return tr[end], W
end

using ProgressMeter

d = 5.0

a_min = 0.0
a_max = 1.0
na = 100
a_vec = LinRange(a_min, a_max, na)

ω_min = 0.2
ω_max = 2.0
nω =  1000
ω_vec = LinRange(ω_min, ω_max, nω)

ds = Systems.vanderpol(; μ = d)
  
wmat = zeros(nω,na)

Ttr =  5000 # transient time
T =   4000 #  time on the attractor used for estimating winding no.
u0 = ds.u0


W_vec = zeros(nω)
@showprogress 1 "Computing..." for ia = 1:na
    a = a_vec[ia]
    for iω = 1:length(ω_vec)
        ω = ω_vec[iω]
        u0, W = winding_number(ds, ω, a, u0, Ttr, T)
        wmat[iω,ia] = 1/W
    end # of frequency loop
end  # of amplitude loop

# %% Plot with all numbers
fig = figure(figsize=(0.6*figx,figy)) 
im1 = plt.pcolormesh(ω_vec,a_vec,wmat',cmap = "gnuplot" )

xlabel(L"\omega") 
ylabel(L"a") 

ax = gca()
ax.set_xlim([ω_min,ω_max])

cbar = plt.colorbar(im1)  # , ticks=[-0.2, 0.2, 0.6, 1.], orientation="horizontal")
cbar.set_label(L"W") 

fig.tight_layout(pad=0.3)

# %% Plot version with only wanted numbers
figure()
qmat = 4*ones(nω,na)   # array with winding numbers to be displayed

rat_num = [1/2, 2/3,   1, 4/3, 3/2, 5/3, 2,  7/3, 5/2, 8/3,  3,  10/3,  7/2]
nrn = length(rat_num)
for ia = 1:na
    for iomega = 1:nω
        wind_no = wmat[iomega,ia] 
        # check whether equal p/q
        for rat in rat_num
            if abs(rat - wind_no) < 0.008  # this threshold should be as small as possible
               qmat[iomega,ia] = rat
            end
        end  
    end
end

im1 = plt.pcolormesh(ω_vec,a_vec,qmat',cmap = "CMRmap") # a color map which ends with white
plt.clim(0.,4.)
ylabel(L"a") 
xlabel(L"\omega"; labelpad = -10) 
cbar = plt.colorbar(im1)
cbar.set_label(L"W") 
fig.tight_layout(pad=0.3)
wsave(plotsdir("9", "vanderPol_tongues"), fig)
