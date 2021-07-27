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

bd = Billiard(Antidot(Float32[0, 0], 1.0f0, false))

# make particles at the rim
js = [(0.2, sind(30)), (1.5, sind(30 + sqrt(2)/10))]

colors = [to_color(COLORS[i]) for i in 1:length(js)]

ps = Particle{Float32}[]
for (ξ, s) in js
    x, y = from_bcoords(Float32(ξ), Float32(s), bd)
    push!(ps, Particle(x..., atan(y[2], y[1])))
end


fig, bmapax = billiard_bmap_plot(bd, ps;
colors = colors, tail = 100000, steps = 100003, backgroundcolor = RGBf0(1,1,1),
ms = 10, vr = 0.1)
GLMakie.ylims!(bmapax, 0.2, 0.9)


# %% Chaotic billiard (sinai)
InteractiveDynamics.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveDynamics.obfill(::Antidot) = RGBAf0(0,0,0,0)
InteractiveDynamics.obls(::Antidot) = nothing

colors = [to_color(COLORS[i]) for i in (5, 1,2)]

bd = billiard_sinai(0.25, 1, 1)

ps = particlebeam(0.2, 0.75, π + π/4 + 0.414235, 100, 0.002)
fig, bmapax = billiard_bmap_plot(bd, ps; colors = colors,
tail = 3100, steps = 3400, backgroundcolor = RGBf0(1,1,1),
ms = 10, vr = 0.05)