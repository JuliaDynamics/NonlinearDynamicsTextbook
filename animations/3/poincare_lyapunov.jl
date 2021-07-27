# %% Interactive GALI psos for henon heiles
# color coded with lyapunov exponent
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using DynamicalSystems, InteractiveDynamics
using OrdinaryDiffEq
import GLMakie

diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

hh = Systems.henonheiles()

potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
E = energy(get_state(hh)...)

function complete(y, py, x)
    V = potential(x, y)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    ic = [x, y, px, py]
    return ic
end

plane = (1, 0.0) # first variable crossing 0

cmap = cgrad(:cividis)
function λcolor(u)
    λ = lyapunov(hh, 10000; u0 = u, diffeq...)
    λmax = 0.06
    v = clamp(λ/λmax, 0, 1)
    return cmap.colors[round(Int, v*255 + 1)]
    # RGBf0(v/2, 0, v)
end

fig, state = interactive_poincaresos(
    hh, plane, (2, 4), complete;
    labels = ("q₂" , "p₂"), color = λcolor, diffeq...
)

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
