using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveDynamics
using DataStructures: CircularBuffer

# Basic figure layouting
fig = Figure(resolution = (1600, 500), figure_padding = 50, dpi = 600)
display(fig)
axs = [fig[i, 2] = Axis(fig) for i in 1:3]
ax_orig = fig[:, 1] = Axis3(fig)
ax_embe = fig[:, 3] = Axis3(fig)
# Link x axis of middle plots
for i in 1:2
    linkxaxes!(axs[3], axs[i])
    hidexdecorations!(axs[i]; grid = false)
end
# Set up a nice projection for Lorenz system
ax_orig.azimuth = 5.855
ax_embe.azimuth = 5.855
ax_orig.elevation = 0.46
ax_embe.elevation = 0.46
ax_orig.zlabeloffset = 20

# Generate a very large timeseries. The animation will run through the timeseries
Δt = 0.01
tail = 800 # in units of Δt
N = 5 # how many tails we will have
τ = 20 # delay time
ds = Systems.lorenz()
traj_orig = trajectory(ds, N*tail*Δt; dt = Δt, Ttr = 1000)
x = traj_orig[:, 1]
t = 0:Δt:N*tail*Δt # time axis
traj_embed = embed(x, 3, τ)

# Set up observables for original
traj_orig_obs = Observable(CircularBuffer{Point3f0}(tail))
append!(traj_orig_obs[], traj_orig[1:tail])
c0 = to_color(COLORSCHEME[1]) # color for 3D attractors
traj_colors = [RGBAf0(c0.r, c0.g, c0.b, i/tail) for i in 1:tail]
orig_endpoint = Observable(Point3f0(traj_orig[tail]...))

# Setup observables for middle timeseries and for embedded plot
cs = [to_color(COLORS[i+1]) for i in 1:3]
xs = [@view(x[(1 + (i-1)*τ):(end - (3-i)*τ)]) for i in 1:3]
t = t[1:length(xs[1])]
xs_points = [Observable(Point2f0(0.0, 0.0)) for i in 1:3]
traj_embed_obs = Observable(CircularBuffer{Point3f0}(tail))
append!(traj_embed_obs[], traj_embed[1:tail])
embed_endpoint = Observable(Point3f0(0,0,0))


# Plot initial trajectory
proj_line = lines!(ax_orig, xprojection; color = cs[1], linewidth = 1)
proj_line.visible = false
lines!(ax_orig, traj_orig_obs; color = traj_colors, linewidth = 2, transparency=true)
scatter!(ax_orig, orig_endpoint; markersize = 4000, color = c0)
# PLot delayed timeseries (animation will simply change time axis range)
for i in 1:3
    lines!(axs[i], t, xs[i]; color = cs[i], linewidth = 2)
    scatter!(axs[i], xs_points[i]; color = cs[i])
end
xlims!(axs[3], (0, tail*Δt))
# Plot reconstructed attractor
lines!(ax_embe, traj_embed_obs; color = traj_colors, linewidth = 2, transparency=true)
scatter!(ax_embe, embed_endpoint; markersize = 4000, color = c0)

# Update/animation function
function update_to_time!(j)
    append!(traj_orig_obs[], traj_orig[j:(tail+j-1)])
    append!(traj_embed_obs[], traj_embed[j:(tail+j-1)])
    # these two lines are needed to actually trigger update.
    traj_orig_obs[] = traj_orig_obs[]
    traj_embed_obs[] = traj_embed_obs[]
    orig_endpoint[] = traj_orig_obs[][end]
    embed_endpoint[] = traj_embed_obs[][end]
    # update points of timeseries
    for i in 1:3
        offset = tail # - (i-1)*τ
        xs_points[i][] = Point2f0((j-1 + offset)*Δt, xs[i][j + offset])
    end
    xlims!(axs[3], ((j-1 + tail÷2)*Δt, (j-1+tail + tail÷2)*Δt))
end
update_to_time!(1)

# Final labelling
ax_embe.xlabel = "x(t)"
ax_embe.ylabel = "x(t + τ)"
ax_embe.zlabel = "x(t + 2τ)"
axs[1].ylabel = "x(t)"
axs[2].ylabel = "x(t + τ)"
axs[3].ylabel = "x(t + 2τ)"
ax_orig.title = "original set"
axs[1].title = "measured timeseries"
ax_embe.title = "reconstructed set"
for prop in ("tick", "label", "ticklabel")
    for (i, f) in enumerate('x':'z')
        @eval setproperty!(ax_embe, Symbol($f, $prop, :color), cs[$i])
    end
end

# %% Perform animation
# for i in 2:(N-2)*tail
#     update_to_time!(i)
#     yield()
#     # sleep(0.001)
# end

record(fig, projectdir("animations", "6", "embedding.mp4"), 1:(N-2)*tail;
    framerate = 60) do i
    update_to_time!(i)
end
