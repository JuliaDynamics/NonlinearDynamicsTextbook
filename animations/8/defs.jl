using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalBilliards, InteractiveDynamics

# Set style to book colors
InteractiveDynamics.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveDynamics.obcolor(::Obstacle) = to_color(COLORS[1])
InteractiveDynamics.obfill(o::Antidot) = RGBAf0(0,0,0,0)
InteractiveDynamics.obls(::Antidot) = nothing
