# %% Delay time impact
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveDynamics

figure = Figure(resolution = (1000, 800))
display(figure)
ax = figure[1, :] = Axis3(figure)
sll = labelslider!(figure, "τ =", 1:100)
figure[2, :] = sll.layout
τ = sll.slider.value

ds = Systems.lorenz()
x = trajectory(ds, 100; Ttr = 100)[:, 1]

R = embed(x, 3, τ[])
Robs = Observable(R.data)

js = (206, 946)

Pobs = Observable([Point3f0(Robs[][j]) for j in js])

lines!(ax, Robs; color = COLORS[1])
scatter!(ax, Pobs; color = :red, markersize = 1000)
ax.azimuth = 5.555530633326979
ax.xlabel = "x(t)"
ax.ylabel = "x(t+τ)"
ax.zlabel = "x(t+2τ)"

on(τ) do t
    Robs[] = embed(x, 3, t).data
    Pobs[] = [Point3f0(Robs[][j]) for j in js]
end

# record_interaction(figure, string(@__FILE__)[1:end-2]*"mp4"; total_time = 15)
