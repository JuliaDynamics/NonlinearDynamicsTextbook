# fitzhugh nagumo figure
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

fig, axs = subplots(2,3; figsize = (figx, 2figy))
fig.tight_layout()

a = 3.0
b = 0.2
ε = 0.01
I = 0.0

ds = Systems.fitzhugh_nagumo([0,0.0]; a, b, ε, I)

xgrid = -0.3:0.02:1.1
ygrid = -0.1:0.01:0.4

function add_nullclines!(ax, a, b, ε, I)
    u = -0.3:0.01:1.1
    w1 = @. u
    w2 = @. a*u*(u-b)*(1-u) + I
    ax.plot(u, w1, "--", color = "C3", lw = 2.0)
    ax.plot(u, w2, ":", color = "C5")
    ax.set_ylim(ygrid[1], ygrid[end])    # why is this necessary?
end

function add_streamlines!(ax, a, b, ε, I)
    ux = zeros(length(xgrid), length(ygrid))
    uy = copy(ux)
    set_parameter!(ds, [a, b, ε, I])
    for (i, x) in enumerate(xgrid)
        for (j, y) in enumerate(ygrid)
            ux[i, j], uy[i, j] = ds.f(SVector(x, y), ds.p, 0)
        end
    end
    ax.streamplot(Vector(xgrid), Vector(ygrid), ux', uy';
        linewidth = 0.5, density = 0.5, color = "C4"
    )
    ax.set_ylim(ygrid[1], ygrid[end])
end

# Fig A (top left)
# ================
ax = axs[1]
add_nullclines!(ax, a, b, ε, I)
add_streamlines!(ax, a, b, ε, I)

# Add tex to axs [1] explaining nullclines
sss = 26
props = Dict(:boxstyle=>"round", :alpha=>1., :facecolor => "white")
ax.text(0.5, 0.33, "\$\\dot{w}>0\$", color = "C3", size = sss, bbox=props)
ax.text(-0.1, 0.33, "\$\\dot{w}<0\$", color = "C3", size = sss, bbox=props)
ax.text(-0.05, 0.1, "\$\\dot{u}<0\$", color = "C5", size = sss, bbox=props)
ax.text(0.45, -0.05, "\$\\dot{u}>0\$", color = "C5", size = sss, bbox=props)
ax.set_ylabel("\$w\$"; labelpad = -20)
ax.set_xlabel("\$u\$"; labelpad = -20)
ax.set_xticks([0, 1])
ax.set_yticks([0, 0.4])

# Fig B & C: two trajectries & zoomin
# ==================
for ax in (axs[3], axs[5])
    add_nullclines!(ax, a, b, ε, I)

    u1 = [0.19, 0.]
    tr1 = trajectory(ds, 200.0, u1)
    ax.plot(columns(tr1)...; color = "C1", lw = 2)
    ax.scatter(tr1[1]...; color = "C1", s = 20)

    u2 = [0.21, 0.]
    tr2 = trajectory(ds, 300.0, u2)
    ax.plot(columns(tr2)...; color = "C0", lw = 2)
    ax.scatter(tr2[1]...; color = "C0", s = 20)

    # ic = (tr1[1] + tr2[1]) ./ 2
    # ax.scatter(ic...; color = "k", s = 10)

    ax.plot(0, 0; marker = "o", mec = "C2", mew = 1,
        markersize = 10, mfc = :white, zorder = 99)
    # ax.axhline(0; color = "k", lw = 0.5, zorder = 1)
end
ax = axs[3]
ax.set_ylabel("\$w\$"; labelpad = -20)
ax.set_xlabel("\$u\$"; labelpad = -20)
ax.set_xticks([0, 1])
ax.set_yticks([0, 0.4])

zbox = ((-0.05, -0.04), (0.25, 0.04))
axis_zoomin!(axs[5],axs[3],  zbox, zbox, "C2"; lw = 1.0)
axs[5].set_xlim(zbox[1][1], zbox[2][1])
axs[5].set_ylim(zbox[1][2], zbox[2][2])
axs[5].set_xticks([])
axs[5].set_yticks([])

# Fig D: timeseries, pulses, refractory period
using OrdinaryDiffEq
pulses_start = [20, 80, 170] # callback times
pulse_width = 4
pulses_end = pulses_start .+ pulse_width
pulses = sort!(vcat(pulses_start, pulses_end))
Ipulse = 0.2
condition(u,t,integ) = t ∈ pulses # pulse times
function affect!(integ)
    i = integ.t ∈ pulses_start ? Ipulse : 0.0
    integ.p[4] = i
end
cb = DiscreteCallback(condition, affect!)

Tf = 250.0
dt = 0.1
prob = ODEProblem(ds, (0.0, Tf))
sol = solve(prob, Tsit5(); callback=cb, tstops = pulses, dtmax = 0.1)

axs[2].plot(sol.t, sol[1, :])
pt = [any(x -> x ≤ t ≤ x + pulse_width, pulses_start) ? Ipulse : 0.0 for t in sol.t]
axs[2].plot(sol.t, pt; color = "C2")
axs[2].set_xlabel("\$t\$"; labelpad = -20)
axs[2].set_ylabel("\$u\$"; labelpad = -20)
axs[2].set_yticks([0, 1])
axs[2].set_xlim(0, Tf)
axs[2].set_xticks([pulses_start..., Tf])

# Fig E: multi stable
a = 8.
b = 0.2
ε = 0.01
I = 0.0

xgrid = -0.3:0.02:1.1
ygrid = -0.2:0.01:1

add_nullclines!(axs[4], a, b, ε, I)
set_parameter!(ds, [a,b,ε,I])
ax = axs[4]

u1 = [0.3, 0.]
tr1 = trajectory(ds, 300.0, u1)
ax.plot(columns(tr1)...; color = "C0")
ax.scatter(tr1[1]...; color = "C0", s = 20)

u2 = [0.55, 0.8]
tr2 = trajectory(ds, 300.0, u2)
ax.plot(columns(tr2)...; color = "C1")
ax.scatter(tr2[1]...; color = "C1", s = 20)

ax.set_ylabel("\$w\$"; labelpad = -20)
ax.set_xlabel("\$u\$"; labelpad = -20)
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])

# Fig F: limit cycle
# -------------------------
ax = axs[6]
a = 3.
b = -0.05
ε = 0.01
I = 0.
set_parameter!(ds, [a,b,ε,I])
add_nullclines!(ax, a, b, ε, I)

tr = trajectory(ds, 400.0, [0.0, 0.1])
ax.plot(columns(tr)...; color = "C0")
ax.scatter(tr[1]...; color = "C0", s = 20)
ax.set_xlim(-0.5,1.1)
ax.set_ylim(-0.1,0.6)
ax.set_ylabel("\$w\$"; labelpad = -20)
ax.set_xlabel("\$u\$"; labelpad = -20)
ax.set_xticks([0, 1])
ax.set_yticks([0, 0.6])

add_identifiers!(fig)
fig.subplots_adjust(bottom = 0.08, left = 0.05, top = 0.95, right = 0.97, hspace = 0.2)
wsave(plotsdir("2", "fitzhugh"), fig)
