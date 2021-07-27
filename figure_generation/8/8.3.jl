using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using InteractiveDynamics, Random
using DynamicalBilliards
import GLMakie
using GLMakie: to_color, RGBf0, RGBAf0
using DynamicalSystems

# Set style to book colors
InteractiveDynamics.obcolor(::Antidot) = to_color(COLORS[3])
InteractiveDynamics.obcolor(::Obstacle) = to_color(COLORS[3])
InteractiveDynamics.obfill(::Antidot) = RGBAf0(0,0,0,0)
InteractiveDynamics.obls(::Antidot) = nothing

w = 0.2f0
bd = billiard_mushroom(1f0, w, 1f0, 0f0)

pc = MushroomTools.randomchaotic(1f0, w, 1f0)
pc = Particle(-0.01f0, 0.2f0, sqrt(3f0)/2)

pr = Particle(0.0f0, 1.2f0, 0.0f0)
pr2 = Particle(0.0f0, 1.9f0, 0.0f0)

ps = [pc, pr, pr2]
colors = [to_color(COLORS[i]) for i in (1,2,4)]

fig, bmapax = billiard_bmap_plot(bd, ps; res = (1600, 600),
colors = colors, tail = 100000, steps = 1000003, backgroundcolor = RGBf0(1,1,1),
ms = 5)

# interactive_billiard_bmap(bd)

# bmapax.xticklabelsize = 32
# bmapax.yticklabelsize = 32
bmapax.xlabelsize = 40
bmapax.ylabelsize = 40
bmapax.ylabelpadding = 15
# GLMakie.save(plotsdir("mushroom.png"), fig)
