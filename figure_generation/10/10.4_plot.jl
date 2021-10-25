# MSF_ploty_v3.jl
# ==================

# master stability functions for the Roessler system with different coupling terms

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
#using DynamicalSystems 
using PyPlot, OrdinaryDiffEq
using StaticArrays
using LinearAlgebra

# coupling matrix
# ---------------

wmat =  [0   7    0    5    0  ;
         7   0    3    0    4  ;
         0   3    0    0    6  ;
         5   0    0    0    0  ;
         0   4    6    0    0 ]

gmat = - wmat

dg = sum(wmat,dims=2)

for i = 1:5
    gmat[i,i] = dg[i]
end

#check row sum
# println( sum(gmat, dims=2) )

# eigenvalues 
evvec = eigvals(gmat)
# println(evvec)



# load results
# -------------
filepath = datadir("master_stability_function.jld")

alpha_vec = load(filepath,"alpha_vec")

LE_vec_xx = load(filepath,"LE_vec_xx")
LE_vec_xy = load(filepath,"LE_vec_xy")
LE_vec_xz = load(filepath,"LE_vec_xz")

LE_vec_yx = load(filepath,"LE_vec_yx")
LE_vec_yy = load(filepath,"LE_vec_yy")
LE_vec_yz = load(filepath,"LE_vec_yz")

LE_vec_zx = load(filepath,"LE_vec_zx")
LE_vec_zy = load(filepath,"LE_vec_zy")
LE_vec_zz = load(filepath,"LE_vec_zz")


# plot version selection
# ======================
fig = figure("MSF selection 1",figsize=(0.5*figx,1.1*figy)) 


# coupling signal x

# Jh = [1 0 0 ; 0 0 0 ; 0 0 0 ]   # x drives x-ODE
plot(alpha_vec,LE_vec_xx, color = COLORS[1], linestyle="-") 

# Jh = [0 0 0 ; 1 0 0 ; 0 0 0 ]  # x drives y-ODE
# plot(alpha_vec,LE_vec_xy, color = COLORS[1], linestyle=":") 

# Jh = [0 0 0 ; 0 0 0 ; 1 0 0 ]  # x drives z-ODE
plot(alpha_vec,LE_vec_xz, color = COLORS[2], linestyle="--") 

# coupling signal y

# Jh = [0 1 0 ; 0 0 0 ; 0 0 0 ]   # y drives x-ODE
# plot(alpha_vec,LE_vec_yx, color = COLORS[4]) 

# Jh = [0 0 0 ; 0 1 0 ; 0 0 0 ]  # y drives y-ODE
plot(alpha_vec,LE_vec_yy, color = COLORS[4], linestyle=":") 

# Jh = [0 0 0 ; 0 0 0 ; 0 1 0 ]  # y drives z-ODE
# plot(alpha_vec,LE_vec_yz, color = COLORS[4], linestyle="--") 

# coupling signal z

# Jh = [0 0 1 ; 0 0 0 ; 0 0 0 ] # z drives x-ODE
# plot(alpha_vec,LE_vec_zx, color = COLORS[2]) 

# Jh = [0 0 0 ; 0 0 1 ; 0 0 0 ]  # z drives y-ODE
# plot(alpha_vec,LE_vec_zy, color = COLORS[2], linestyle=":") 

# Jh = [0 0 0 ; 0 0 0 ; 0 0 1 ]  # z drives z-ODE
plot(alpha_vec,LE_vec_zz, color = COLORS[3], linestyle="--") 



xlabel(L"$\alpha$";labelpad = -20)
ylabel(L"$\lambda_1$") 

LE_min = -0.3
LE_max = 0.3

ax = gca()
ax.set_xlim([alpha_vec[1], alpha_vec[end]])
ax.set_ylim([LE_min,LE_max])
ax.legend([L"x \to \dot x"; L"x \to \dot z"; L"y \to \dot y"; L"z \to \dot z" ], fontsize = 28,
           loc = "upper right", ncol = 2, handlelength = 1. , handletextpad = 0.5)


# zero line
plot(alpha_vec, 0 .* LE_vec_xx, linewidth=1.0,linestyle="-",color="k" )
ylabel(L"$\lambda_{1}$") #

# plot eigenvalues
sigma = 0.15
scatter(sigma*evvec, 0 .* evvec,color="k",s = 40. )

#add_identifiers!(fig)
fig.tight_layout(pad=0.3)


wsave(plotsdir("10", "MSF_Roessler_selection"), fig)