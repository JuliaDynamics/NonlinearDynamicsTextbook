include("defs.jl")

bd = Billiard(Antidot([0.0, 0.0], 1.0, false))

js = [(0.2, sind(30)), (1.5, sind(30 + sqrt(2)/10))]
colors = [to_color(COLORS[i+1]) for i in 1:length(js)]

ps = Particle{Float64}[]
for (ξ, s) in js
    x, y = from_bcoords(ξ, s, bd)
    push!(ps, Particle(x..., atan(y[2], y[1])))
end

fig, = interactive_billiard(
    bd, ps;
    colors, tail = 20000, backgroundcolor = RGBf0(1,1,1),
)


record_interaction(fig, joinpath(@__DIR__, "circle_billiard.mp4"))
