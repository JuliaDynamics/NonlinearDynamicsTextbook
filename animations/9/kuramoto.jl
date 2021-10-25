# %% Interactive bifurcations for 1D energy balance
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using DynamicalSystems, OrdinaryDiffEq, Random
using InteractiveDynamics: record_interaction
import GLMakie

D = 25
omega_vec = randn(D)

ds = Systems.kuramoto(D)
diffeq = (alg = Tsit5(), adaptive = false, dt = 0.05)
integ = integrator(ds; diffeq...)

fig = GLMakie.Figure()
display(fig)
cmap = :curl
axku = GLMakie.Axis(fig[1,1])


# Other static elements
GLMakie.hidedecorations!(axku)
axku.aspect = GLMakie.DataAspect()
GLMakie.lines!(axku, cos.(0:0.001:2π+0.002), sin.(0:0.001:2π+0.002);color = :black)

# Plot balls
phases = GLMakie.Observable(copy(integ.u))
balls = GLMakie.lift(u -> [GLMakie.Point2f0(cos(φ), sin(φ)) for φ in u], phases)
GLMakie.scatter!(axku, balls;
    markersize  = 12, color = 1:D, colormap = cmap, strokewidth=2, strokecolor = :black
)

# Plot arrows
arrows_end = GLMakie.lift(b -> 0.8b, balls)
arrows_start = [GLMakie.Point2f0(0,0) for i in 1:D]
GLMakie.arrows!(axku, arrows_start, arrows_end; 
    color = 1:D, colormap = (cmap, 0.5), arrowsize = 10, linewidth = 4,
)

# Plot mean field Vector
# add sliders and buttons
fig[2, 1] = controllayout = GLMakie.GridLayout(tellwidth = false)
run = controllayout[1, 1] = GLMakie.Button(fig; label = "run")
Ks = 0.0:0.1:10.0
kslider = GLMakie.labelslider!(fig, "K =", Ks; sliderkw = Dict(:startvalue => ds.p.K))
controllayout[1, 2] = kslider.layout


# run button functionality
isrunning = GLMakie.Observable(false)
GLMakie.on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
GLMakie.on(run.clicks) do clicks
    @async while isrunning[]
        GLMakie.isopen(fig.scene) || break # ensures computations stop if closed window
        step!(integ)
        phases[] = integ.u
        yield()
        sleep(0.001) # or `yield()` instead
    end
end

GLMakie.on(kslider.slider.value) do val
    integ.p.K = val
    u_modified!(integ, true)
end

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 15, sleep_time = 2)
