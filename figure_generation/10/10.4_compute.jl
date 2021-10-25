# MSF_compute_v3.jl
# ==================

# master stability functions for the Roessler system with different coupling terms

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
#using DynamicalSystems 
using PyPlot, OrdinaryDiffEq
using StaticArrays


function Roessler_ODE(u,p,t)
    a, b, c = p
    x, y, z = u
    dx =  - y - z
    dy =  x + a*y
    dz = b + z*(x-c)
    return SVector(dx, dy, dz)
end


function Roessler_MSF_ODE(u,p,t)
    a, b, c, alpha, Jh = p
    x1, x2, x3, y1, y2, y3  = u
    dx1 =  - x2 - x3
    dx2 =  x1 + a*x2
    dx3 = b + x3*(x1-c)

   # linearized equation
    Jg = [0 -1 -1 ; 1 a 0 ; x3 0 x1-c]
    dy = (Jg .- alpha.*Jh) * u[4:6]

    return SVector(dx1, dx2, dx3, dy[1], dy[2], dy[3])
end

function MSF(uinit,a,b,c,Jh)
  
    nle = 20000 # number of renormalization steps
    tle = 6.  # period of time between renormalizations

    alpha_min = -2.
    alpha_max = 8.
    nalpha = 200
    alpha_vec = LinRange(alpha_min, alpha_max, nalpha)
  #  avec = zeros(nalpha)   # alpha values
    LE_vec = zeros(nalpha)  # corresponding Lyapunov exponent
    
   # uinit = u0[1:3]
    
    for ialpha = 1:nalpha
    
        alpha = alpha_vec[ialpha] 
        cum = 0.

        p = a, b, c, alpha, Jh     
        u0 = [uinit ; SVector(1., 0., 0.) ]           

        for ile = 1:nle   # loop over renormalization steps

            tspan = (0., tle)
            prob = ODEProblem(Roessler_MSF_ODE, u0, tspan, p)
            sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,tsamp = tle)
    
            yvec = [u[4:6] for u in sol.u]  # perturbation vector
            yend = yvec[end]
            ynorm = sqrt( sum(yend .* yend)  )
            yend = yend / ynorm      # renormalization
            cum = cum + log(ynorm)  
    
            uvec = [u[1:3] for u in sol.u]
            uend = uvec[end]
    
            u0 = [uend ; yend]   # initial condition for next integration step
    
        end
    
        lambda_max = cum / (nle*tle) # Lyapunov exponent
       # println(alpha, "  ", lambda_max)
    
       # avec[ialpha] = alpha
        LE_vec[ialpha] = lambda_max
    
    end # of alpha loop
    return alpha_vec, LE_vec
    end # of function
    
# ------------------------------------------------------------------------

# main program
# ------------
a = 0.2
b = 0.2 
c = 5.7 
#Jh = [0 0 0 ; 0 1 0 ; 0 0 0 ]


# transient
# ---------
p = a, b, c
u0 = SVector(0.1, 0.1, 0.)

tspan = (0., 1000)
prob = ODEProblem(Roessler_ODE, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,tsamp = 100.)


# attractor  (just to check that the solution is o.k.)
# ---------
u = sol.u
# u0 = [u0[end] ; SVector(1., 0., 0.) ]
uinit = u[end]
println(uinit)

p = a, b, c   #, alpha, Jh

tattr = 200. 
tspan = (0., tattr)
prob = ODEProblem(Roessler_ODE, uinit, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,tsamp = 0.1)

tvec = sol.t
x1vec = [u[1] for u in sol.u]
x2vec = [u[2] for u in sol.u]
x3vec = [u[3] for u in sol.u]

# plot
# ----
fig = figure("Roessler",figsize=(0.5*figx,figy)) 
plot(x1vec,x2vec) 
xlabel(L"$x_1$") 
ylabel(L"$x_2$") 

# plot
# ----
fig = figure("Roessler ts",figsize=(0.5*figx,figy)) 
plot(tvec[1:200],x1vec[1:200]) 
xlabel(L"$t$") 
ylabel(L"$x_1$") 

fig.tight_layout(pad=0.3)


# compute MSF
# -----------
fig = figure("MSF 2",figsize=(figx,figy)) 

LE_min = -0.015
LE_max = 0.02

# coupling signal x
subplot(1,3,1)

Jh = [1 0 0 ; 0 0 0 ; 0 0 0 ]   # x drives x-ODE
alpha_vec, LE_vec_xx = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_xx, color = COLORS[1], linestyle="-") 

Jh = [0 0 0 ; 1 0 0 ; 0 0 0 ]  # x drives y-ODE
alpha_vec, LE_vec_xy = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_xy, color = COLORS[1], linestyle=":") 

Jh = [0 0 0 ; 0 0 0 ; 1 0 0 ]  # x drives z-ODE
alpha_vec, LE_vec_xz = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_xz, color = COLORS[1], linestyle="--") 


# zero line
plot(alpha_vec, 0 .* LE_vec_xx, linewidth=1.0,linestyle="-",color="k" )


xlabel(L"$\alpha$") 
ylabel(L"$\lambda_{max}$") 

ax = gca()
ax.set_xlim([alpha_vec[1], alpha_vec[end]])
ax.set_ylim([LE_min,LE_max])
ax.legend([L"x \to x"; L"x \to y"; L"x \to z"], fontsize = 28, loc = "upper left")


# coupling signal y
subplot(1,3,2)
Jh = [0 1 0 ; 0 0 0 ; 0 0 0 ]   # y drives x-ODE
alpha_vec, LE_vec_yx = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_yx, color = COLORS[4]) 

Jh = [0 0 0 ; 0 1 0 ; 0 0 0 ]  # y drives y-ODE
alpha_vec, LE_vec_yy = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_yy, color = COLORS[4], linestyle=":") 

Jh = [0 0 0 ; 0 0 0 ; 0 1 0 ]  # y drives z-ODE
alpha_vec, LE_vec_yz = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_yz, color = COLORS[4], linestyle="--") 

# zero line
plot(alpha_vec, 0 .* LE_vec_yx, linewidth=1.0,linestyle="-",color="k" )

xlabel(L"$\alpha$") 

ax = gca()
ax.set_xlim([alpha_vec[1], alpha_vec[end]])
ax.set_ylim([LE_min,LE_max])
ax.legend([L"y \to x"; L"y \to y"; L"y \to z"], fontsize = 28)




# coupling signal z
subplot(1,3,3)
Jh = [0 0 1 ; 0 0 0 ; 0 0 0 ] # z drives x-ODE
alpha_vec, LE_vec_zx = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_zx, color = COLORS[2]) 

Jh = [0 0 0 ; 0 0 1 ; 0 0 0 ]  # z drives y-ODE
alpha_vec, LE_vec_zy = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_zy, color = COLORS[2], linestyle=":") 

Jh = [0 0 0 ; 0 0 0 ; 0 0 1 ]  # z drives z-ODE
alpha_vec, LE_vec_zz = MSF(uinit,a,b,c,Jh)
plot(alpha_vec,LE_vec_zz, color = COLORS[2], linestyle="--") 

# zero line
plot(alpha_vec, 0 .* LE_vec_zx, linewidth=1.0,linestyle="-",color="k" )

xlabel(L"$\alpha$") 

ax = gca()
ax.set_xlim([alpha_vec[1], alpha_vec[end]])
ax.set_ylim([LE_min,LE_max])
ax.legend([L"z \to x"; L"z \to y"; L"z \to z"], fontsize = 28)


add_identifiers!(fig)
fig.tight_layout(pad=0.3)

# save
# ----
# dname = "/Users/parlitz/ownCloud/NonlinearDynamicsTextbook/figure_generation/10/MSF/"
# fname = "MSF_Roessler_6_20000"
# ftype = ".png"
# savefig(dname*fname*ftype)

save(datadir("master_stability_function.jld"),"alpha_vec",alpha_vec,
"LE_vec_xx",LE_vec_xx, "LE_vec_xy",LE_vec_xy,"LE_vec_xz",LE_vec_xz,
"LE_vec_yx",LE_vec_yx, "LE_vec_yy",LE_vec_yy,"LE_vec_yz",LE_vec_yz,
"LE_vec_zx",LE_vec_zx, "LE_vec_zy",LE_vec_zy,"LE_vec_zz",LE_vec_zz)

