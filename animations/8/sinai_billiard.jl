include("defs.jl")

bd = billiard_sinai(0.25)
ps = particlebeam(0.2, 0.75, π + π/4 + 0.414235, 100, 0.002)

fig, = interactive_billiard(
    bd, ps;
    tail = 500, backgroundcolor = RGBf0(1,1,1),
    vr = 0.05
)

# record_interaction(fig, joinpath(@__DIR__, "sinai_billiard.mp4"))


# %% Ergodicity in the Sinai billiard
bd = billiard_sinai(0.25)

fig, = interactive_billiard_bmap(bd)

# record_interaction(fig, joinpath(@__DIR__, "sinai_ergodicity.mp4"); total_time = 20)
