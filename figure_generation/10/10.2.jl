using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("colorscheme.jl"))
import GLMakie, GraphMakie
using GraphMakie.NetworkLayout, LightGraphs

n0 = 25
seed = 42
n2 = 100 # for getting distribution


er = n -> erdos_renyi(n, 0.25; seed)
ba = n -> barabasi_albert(2n, 3, 1; complete=true, seed)
ws = n -> watts_strogatz(n, 4, 0.1; seed)

fig = GLMakie.Figure(resolution = (800, 400))

graphs = (er, ba, ws)
names = ("Erdős–Rényi ", "Barabási–Albert", "Watts-Strogatz")
for (i, n) in enumerate(names)
    ax = GLMakie.Axis(fig[1,i], title = n)
    g = graphs[i](n0)
    layout = i == 3 ? Shell() : Spring()
    GraphMakie.graphplot!(ax, g, layout=layout, edge_color = COLORS[1], node_color = BLACK)
    GLMakie.hidedecorations!(ax)
    GLMakie.hidespines!(ax)

    degax = GLMakie.Axis(fig[2,i], ylabel = i == 1 ? L"P(k)" : "")
    degax.ylabelsize = 32
    g2 = graphs[i](n2)
    dh = degree_histogram(g2)
    k = Int.(keys(dh))
    p = Int.(values(dh))
    GLMakie.barplot!(degax, k, p, color = BLACK)
end

GLMakie.rowsize!(fig.layout, 2, GLMakie.Relative(1/3))
# display(fig)
GLMakie.save(plotsdir("10", "network_showcase.png"), fig; px_per_unit = 4)