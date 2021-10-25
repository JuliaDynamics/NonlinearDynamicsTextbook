using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# sensitive dependence
Random.seed!(5)
fig = figure(figsize=(figx/2, figy))
ax = gca()
u0 = [10,10.0,10]
lo = Systems.lorenz([10,10.0,10])

for i in 1:3
    u = u0 .+ i*1e-3
    tr = trajectory(lo, 15, u)
    plot(0:0.01:15, tr[:, 1]; c = COLORS[[1,4,3][i]])
end

xlabel("\$t\$",labelpad = -15)
yticks(-15:10:15)
xticks(0:5:15)
ylabel("\$x\$",labelpad = -20)
tight_layout(pad = 0.25)
wsave(plotsdir("3","sensitive"), fig)
