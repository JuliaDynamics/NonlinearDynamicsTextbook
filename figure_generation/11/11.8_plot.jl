# Aliev_Panfilov_2D_plot_sequence_v1.jl
# =====================================   u.p. 5.8.21  / 5.9.21

# Aliev-Panfilov model
# load data and plot sequence with phase singularities

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using StaticArrays
using JLD
using Random
using Statistics
using PyPlot

# load data
# ---------    
path = datadir("Aliev_Panfilov_snapshots.jld")
L = load(path,"L")
N = load(path,"N")
tvec = load(path,"tvec")
uout = load(path,"uout")   

h = L/(N-1) # spatial steps size
xx = h * (0:N-1)   # x coordinates of grid
yy = h * (0:N-1)   # y coordinates of grid
 
nts = length(tvec)
println("number of snapshots nts = ", nts)

# plot selection

nsel =  [16, 18, 20, 22]  # 75, 85, 95, 105
N = length(nsel)
# --------------
fig, axs = subplots(1, N; figsize = (figx, 0.8figy))

# nsel =  [30, 32, 34, 36]  # 145, 155, 165, 175
cmap = "BuPu"
for iplt = 1:N
  
    n = nsel[iplt]
    ures = uout[n]

    uu = ures[1,:,:]
    vv = ures[2,:,:] 

    subplot(1,4,iplt)
    im = plt.pcolormesh(xx,yy,uu; cmap)  # "gnuplot" ) #"inferno")

    ax = gca()
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    ax.set_aspect("equal")

   # plt.title(L"t = "*string(floor(Int,tvec[n])), fontsize = 24)
    println("t=", tvec[n])
    # if iplt == N
    #     # cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.9])
    #     # cb = colorbar(im; cax = cbar_ax)
    #     fig.colorbar(im, ax=axs.ravel().tolist())
        
    # end
end

fig.tight_layout(pad=0.3)
# fig.subplots_adjust(right=0.92)
wsave(plotsdir("11", "Aliev_Panfilov_sequence"), fig)
