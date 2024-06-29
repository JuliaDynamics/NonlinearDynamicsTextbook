using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("theme.jl"))
using CairoMakie
using StaticArrays

xmax = 1.0
xmin = -1.0
ymin = -1.0
ymax = 1.0
d = 50
const a = 1.4
const b = 0.3

ics = [
    SVector(x, y) for x in range(xmin, xmax; length = d) for y in range(ymin, ymax; length = d)
]

f1(u) = SVector(u[1], 1 - a*u[1]^2 + u[2])
f2(u) = SVector(b*u[1], u[2])
f3(u) = SVector(u[2], u[1])
function three_phases(ics)
    p1 = f1.(ics)
    p2 = f2.(p1)
    p3 = f3.(p2)
    return p1, p2, p3
end

function linear_transition!(o, p1, p2, steps, io)
    for i in 1:steps # assumes we start with o[] = p1
        o[] = p1 .+ (p2 .- p1) .* i/steps
        recordframe!(io)
    end
    recordframe!(io)
end

iterations = 10
steps = 20

fig = Figure(); display(fig)
ax = Axis(fig[1,1])
o = Observable(ics)
scatter!(ax, o;
    strokewidth = 0.5, strokecolor = :black, markersize = 5
)
ax.limits = ((-2.5,2.5), (-1.5,1.5))

record(fig, string(@__FILE__)[1:end-2]*"mp4"; framerate = 30) do io
    for j in 1:iterations
        p1, p2, p3 = three_phases(ics)
        linear_transition!(o, ics, p1, steps, io)
        linear_transition!(o, p1, p2, steps, io)
        linear_transition!(o, p2, p3, steps, io)
        global ics = p3
        recordframe!(io)
    end
end