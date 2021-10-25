using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))

using InteractiveDynamics
using DynamicalSystems, OrdinaryDiffEq
import GLMakie

function rossler_driving(u, p, t)
    k = p[1]
    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3 = u

    dx_1  =  - x_2 -x_3 
    dx_2  =   x_1 + 0.2x_2  
    dx_3  =  0.2 + x_3*(x_1-5.7)
    (dx_1, dx_2, dx_3) = 6 .* (dx_1, dx_2, dx_3)
    
    dy_1  =  10*(- y_1 + y_2)
    dy_2  =  28y_1 - y_2 - y_1*y_3 + k*x_2 
    dy_3  =  y_1*y_2 - 2.666y_3
        
    dz_1  =  10*(- z_1 + z_2)
    dz_2  =  28z_1 - z_2 - z_1*z_3 + k*x_2 
    dz_3  =  z_1*z_2 - 2.666z_3

    return SVector(dx_1, dx_2, dx_3, dy_1, dy_2, dy_3, dz_1, dz_2, dz_3)
end

x0 = [-1, -2, 0.1]
y0 = [5.0, 5.0, 20.0]
z0 = [-5.0, -5.0, 20.0]

u0 = vcat(x0, y0, z0)
k = 5.0

ds = ContinuousDynamicalSystem(rossler_driving, u0, [k])

ps = Dict(
    1 => 0.0:0.01:12,
)
pnames = Dict(
    1 => "k",
)

u0s =  [ds.u0]

fig, obs = interactive_evolution_timeseries(
    ds, u0s, ps; tail = 500, idxs = [4,7], pnames,
    lims = ((-20,20), (-20,20))
)

GLMakie.Label(fig[:, :, GLMakie.Top()], "Roessler subsystem driving two Lorenz subsystems: y, z. Plotted are the first variables of each subsystem. See Eq.(9.9)", valign = :bottom, padding = (0, 0, 10, 0))
GLMakie.content(fig[1,1]).xlabel = "y₁"
GLMakie.content(fig[1,1]).ylabel = "z₁"
GLMakie.content(fig[:, 2][1,1]).ylabel = "y₁"
GLMakie.content(fig[:, 2][2,1]).ylabel = "z₁"
GLMakie.content(fig[:, 2][2,1]).xlabel = "t"


# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 20, sleep_time = 2)
