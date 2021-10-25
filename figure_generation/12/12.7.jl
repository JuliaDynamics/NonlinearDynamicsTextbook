# %% Saltzman iceage models
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

function sm91(state, param, t)
    u, ν, p, r, s, q = param
    x, y, z = state
    dx = -x -y -u*z # - ν*astronomical(t/10) # ensure time is in kyr for interpolation
    dy = -p*z +r*y -s*y^2 - y^3
    dz = -q*(x+z)
    return SVector(dx, dy, dz)
end

# Astronomical forcing:
using DelimitedFiles, Statistics
X = readdlm(projectdir("exercise_data", "16.csv"))
tsol = X[:, 1] # time in kyr
Isol = X[:, 5]
# Transform simulation to last 500 kyr and reverse time axis
# L = length(Isol)
# Isol = Isol[1:L÷2]
# tsol = tsol[1:L÷2]
Isol = standardize(Isol)
tsol .+= abs(tsol[1])
using Interpolations
astronomical = LinearInterpolation(tsol, Isol)
quasiperiodic(t) = cospi(2t/22.1) - cospi(2t/44.1)
# Default parameter values
p = 1.0
q = 2.5
r = 1.3
s = 0.6
u = 0.2
ν = 0.0
t₀ = 10.0 # timescale in kyr
params = [u, ν, p, r, s, q]
u₀ = [-0.2, 1.0, 0.0]

ds = ContinuousDynamicalSystem(sm91, u₀, params)

fig, axs = subplots(1,3)

# Unforced solution
T = 50.0
Δt = 0.01
tr = trajectory(ds, T; Δt)
t = (0:Δt:T) .* t₀
for i in 1:3; 
    x = tr[:, i]
    axs[1].plot(t, x, ls = LINESTYLES[i], color = "C$(i-1)")
    imax = argmax(x)
    label = "\$"*("x", "y", "z")[i]*"\$"
    axs[1].text(t[imax], x[imax], label; color = "C$(i-1)")
end
axs[1].set_xlabel("\$t\$ (kyr, arbtr.)")
axs[1].set_ylim(-1.5, 1.5)
fig.tight_layout(pad=0.3)

# # Forced solution
# set_parameter!(ds, 2, 0.5)
# tr = trajectory(ds, 120; Δt)
# axs[3].plot(tr[:, 2])


# transform tsol into same time units as above
# axs[2].plot(tsol, Isol; c = "C3")
# axs[2].set_xlabel("\$t\$ (kyr)")
# axs[2].text(200, 410, "\$F(t)\$", color = "C3")

# Van der poll version
function vdp(u, p, t)
    x, y = u
    α, β, γ, τ, A = p
    dx = -y - β + γ*astronomical_forcing(A, t)
    # dx = -y - β + γ*sinpi(2t/41)
    dy = -α*(y^3/3 - y - x)
    return SVector(dx/τ, dy/τ)
end

include("12.7_astronomical_forcing.jl")



α = 11.11
β = 0.25
γ = 0.75
τ = 40.09
p = (α, β, γ, τ, A)

u₀ = [0.1, -0.2]
T = 600.0
Δt = 0.1
t0 = 500.0
t = (0:Δt:T) .- t0

F = astronomical_forcing.(Ref(A), t)
axs[2].plot(t, F; color = "C3")
axs[2].text(-400, -5, "\$F(t)\$", color = "C3")
# %%
axs[3].clear()

ds = ContinuousDynamicalSystem(vdp, u₀, p; t0 = -t0)
tr = trajectory(ds, T; Δt)
axs[3].plot(t, standardize(tr[:, 1]))

# Plot vostok
V = readdlm(projectdir("exercise_data", "17.csv"))
Vt = V[:, 1] ./ 1000 # in ky
VT = V[:, 2]
Vt = -reverse(Vt)
VT = reverse(VT)
# Interpolate Vostok data for better plot clarity
using Dierckx
# interpolation object:
spl = Spline1D(Vt, VT; k=3, bc="extrapolate")
Vti = Vt[1]:(10Δt):Vt[end]
VTi = spl(Vti)


axs[3].plot(Vti, standardize(VTi); 
    color = "C2", zorder = 99, alpha = 0.8, ls = "-.", lw = 2,
    # marker = "o", mew = 1, mfc = "none", alpha = 0.2, ms = 5,
)

axs[3].text(38.6, 1.61, "\$x\$"; color = "C0")
axs[3].text(-120, 2.6, "Vostok"; color = "C2")

axs[2].set_xlabel("\$t\$ (kyr)")
axs[3].set_xlabel("\$t\$ (kyr)")

add_identifiers!(fig)
wsave(plotsdir("12", "iceages"), fig)
