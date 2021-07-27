include("defs.jl")

# Zooming into standard map
using DynamicalSystems
figure = Figure(resolution = (800,800))
ax = figure[1,1] = Axis(figure)
deactivate_interaction!(ax, :rectanglezoom)
rect = select_rectangle(ax.scene)
sm = Systems.standardmap(; k = 1.0)
g = 10 # grid density
inrect(u, r) = (r.origin[1] ≤ u[1] ≤ r.origin[1]+r.widths[1]) && (r.origin[2] ≤ u[2] ≤ r.origin[2] + r.widths[2])
integ = integrator(sm)
marker = InteractiveDynamics.MARKER

N = 100 # how many points to plot inside the rect
Nslider = labelslider!(
    figure, "N =", round.(Int, range(100, 100000; length = 100))
)
N = Nslider.slider.value

recompute = Button(figure; label = "recompute")

figure[2,1][1,1] = recompute
figure[2,1][1,2] = Nslider.layout

sco = Observable([Point2f0(sm.u0)])

scatter!(ax, sco; markersize = 2*px, strokewidth = 0, color = :black)

on(rect) do r # new subgrid selected
    xs = range(r.origin[1], r.origin[1]+r.widths[1]; length = g)
    ys = range(r.origin[2], r.origin[2]+r.widths[2]; length = g)
    s = sco[]
    n = N[]
    empty!(s)
    for x in xs
        for y in ys
            DynamicalSystems.reinit!(integ, SVector(x, y))
            push!(s, integ.u)
            i = 0
            for i in 1:n
                DynamicalSystems.step!(integ)
                if inrect(integ.u, r)
                    push!(s, integ.u)
                end
            end
        end
    end
    sco[] = s
    xlims!(ax, xs[1], xs[end])
    ylims!(ax, ys[1], ys[end])
end

on(recompute.clicks) do c
    rect[] = rect[]
end

rect[] = FRect2D([0,0],[2π,2π])
hidexdecorations!(ax)
hideydecorations!(ax)

display(figure)

# record_interaction(figure, joinpath(@__DIR__, "sm_zoom.mp4"); total_time = 20)
