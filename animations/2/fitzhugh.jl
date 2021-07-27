# Interactive clickable FitzHughNagumo state space

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveDynamics

a = 3.
b = -0.05
ε = 0.01
I = 0.0

fh = Systems.fitzhugh_nagumo(zeros(2); a, b, ε, I)

# Layout plot:
fig = Figure(resolution = (1000, 800))
ax = fig[1,1] = Axis(fig)
display(fig)
ax.xlabel = "u"
ax.ylabel = "w"
ax.title = "FitzHugh Nagumo, a=$a, b=$b"

# plot the two nullclines
us = -0.4:0.01:1.2
w1 = us
w2 = @. a*us*(us-b)*(1-us) + I
lines!(ax, us, w1; color = COLORS[1], linewidth = 2, linestyle = :dash)
lines!(ax, us, w2; color = COLORS[2], linewidth = 2, linestyle = :dot)

# Create trajectories on point selection
spoint = select_point(ax.scene)
on(spoint) do pos
    tr = trajectory(fh, 10000, SVector(pos...))
    lines!(ax, columns(tr)...;
        color = InteractiveDynamics.randomcolor(), linewidth = 4.0
    )
end

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 10)
