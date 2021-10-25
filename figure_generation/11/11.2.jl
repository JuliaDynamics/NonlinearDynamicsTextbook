# This file needs the source file to have run with appropriate parameters!!!
# It loads and plots the results!
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using JLD, PyPlot

ds = Float64[4, 3, 2]

a = 9.0
b = 9.8
L = 50.0
N = 250
ic = :random
prefix = "brusselator2D"
suffix = "jld2"
fnames = [savename(prefix, @strdict(a, b, d, L, N, ic), suffix) for d in ds]
fig = figure(figsize=(figx, 1.6*figy))

iplt = 0 
nfile = length(fnames)
println(nfile)

for ifile = 1:nfile
    fname = fnames[ifile]
    println(fname)
    data = load(datadir("Brusselator", fname))
    uout = data["u"]
    tvec = data["t"]
    nts = length(tvec)
    for n = 1:nts
        ures = uout[n]
        iplt += 1 
        ax = subplot(2,4,iplt)
        im = ax.pcolormesh(ures[1,:,:]; cmap = "inferno")
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        plt.title(L"t = "*string(floor(Int,tvec[n])), fontsize = 24)
        ax.set_aspect("equal")
    end
end # of iplt loop

# plt.colorbar(im)   
fig.tight_layout(pad=0.4)
add_identifiers!(fig)
wsave(plotsdir("11", "Brusselator_2D_Turing_patterns"), fig)