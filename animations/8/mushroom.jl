include("defs.jl")

w = 0.2f0
bd = billiard_mushroom(1f0, w, 1f0, 0f0)

pc = MushroomTools.randomchaotic(1f0, w, 1f0)
pc = Particle(-0.01f0, 0.2f0, sqrt(3f0)/2)

pr = Particle(0.0f0, 1.2f0, 0.0f0)
pr2 = Particle(0.0f0, 1.9f0, 0.0f0)

ps = [pc, pr, pr2]
colors = [to_color(COLORS[i+1]) for i in 1:length(ps)]

fig, = interactive_billiard(bd, ps; colors = colors, tail = 10000)

# record_interaction(fig, joinpath(@__DIR__, "mushroom.mp4"))
