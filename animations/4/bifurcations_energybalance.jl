# %% Interactive bifurcations for 1D energy balance
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using DynamicalSystems, InteractiveDynamics
using Roots
import GLMakie

αtan(T) = 0.5 - (0.4/π)*atan((T-263)/2)
dTdt(T, ε = 0.65, α=αtan, s= 1.0) = s*(1 - α(T)) - 1.6e-10 * ε * T^4

Ts = 200:0.5:320.0
εs = 0.2:0.01:1.0

function roots(ε)
    f = T -> dTdt(T, ε)
    r = Roots.find_zeros(f, Ts[1], Ts[end])
    s = [GLMakie.to_color(COLORS[isodd(i) ? 4 : 3]) for i in 1:length(r)]
    return [Point2f0(x, 0) for x in r], s # roots, stability
end

fig = GLMakie.Figure(resolution = (800, 500))

ax = fig[1, 1] = GLMakie.Axis(fig;
    xlabel = GLMakie.L"T", ylabel = GLMakie.L"f(T)",
    ylabelsize = 24, xlabelsize = 24
)
ax.limits = ((Ts[1], Ts[end]), (-0.3, 0.3))
axb = fig[1, 2] = GLMakie.Axis(fig;
    xlabel = GLMakie.L"\epsilon", ylabel = GLMakie.L"T^*",
    ylabelsize = 24, xlabelsize = 24
)
sll = GLMakie.labelslider!(fig, "ϵ =", εs; sliderkw = Dict(:startvalue => 0.65))
fig[2, :] = sll.layout
ε_observable = sll.slider.value
axb.limits = ((εs[1], εs[end]), (220, 300))


# initialize plot
fs = GLMakie.Observable(dTdt.(Ts))
GLMakie.lines!(ax, Ts, fs; color = COLORS[2], linewidth = 4)
GLMakie.hlines!(ax, 0; color = COLORS[1])
r, s = roots(0.65)
rootvals = GLMakie.Observable(r)
rootcols = GLMakie.Observable(s)
GLMakie.scatter!(ax, rootvals; color = rootcols, markersize = 10)

stablefp = GLMakie.Observable(GLMakie.Point2f0[])
unstablefp = GLMakie.Observable(GLMakie.Point2f0[])

GLMakie.scatter!(axb, stablefp; color = COLORS[4])
GLMakie.scatter!(axb, unstablefp; color = COLORS[3])

display(fig)

GLMakie.on(ε_observable) do ε
    fs[] = dTdt.(Ts, ε)
    r, s = roots(ε)
    rootvals[] = r
    rootcols[] = s

    for (i, ρ) in enumerate(r)
        p = GLMakie.Point2f0(ε, ρ[1])
        if isodd(i)
            push!(stablefp[], p)
        else
            push!(unstablefp[], p)
        end
    end
    stablefp[] = stablefp[]
    unstablefp[] = unstablefp[]
end

ε_observable[] = 0.65

# %%
record_interaction(fig, string(@__FILE__)[1:end-2]*"mp4"; total_time = 15, sleep_time = 2)
