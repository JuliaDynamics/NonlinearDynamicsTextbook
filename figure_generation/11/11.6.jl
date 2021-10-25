# requires 11.6_produce.jl to have run first
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using PyPlot
using Statistics

function compute_phases(uu,vv,xx,yy) # also computes phase integral
    umean = mean(uu)
    vmean = mean(vv)
    ustd = std(uu)
    vstd = std(vv)  
    N = size(uu, 1)
    
    phimat = atan.( (vv .- vmean) ./ vstd , (uu .- umean) ./ ustd)     

    # compute path integrals
    kk = [0 1 1 0 0]
    ll = [0 0 1 1 0]

    xps = Float64[]    # coordinates of phase singularities (to be detected)
    yps = Float64[]
    for i=1:N-1
        for j = 1:N-1
            sum = 0 
            for m = 1:4
                diffalpha = phimat[i+kk[m+1],j+ll[m+1]] - phimat[i+kk[m],j+ll[m]]
                diffalpha = mod( diffalpha + pi, 2 .* pi) - pi
                sum = sum + diffalpha
            end
                    
            topcharge = sum / (2. * pi)

            if abs(topcharge) > 0.95   # phase singularity detected
                xpos = 0.5*(xx[i] + xx[i+1])
                ypos = 0.5*(yy[j] + yy[j+1])
                push!(xps, xpos)
                push!(yps, ypos)
                # xps = cat(xps, xpos, dims=1)
                # yps = cat(yps, ypos, dims=1) 
                # println("xps = ",xps)
            end
        end
    end
    return phimat, xps, yps
end

data = wload(datadir("FitzHugh", "spiralwave.jld2"))
@unpack uout, tvec, params = data
@unpack a, b, d, ε, hsq6, N, L = params
h = L/N # spatial steps size
xx = h * (0:N-1)   # x coordinates of grid
yy = h * (0:N-1)   # y coordinates of grid

nt = length(tvec)
ures = uout[nt]

uu = ures[1,:,:]
vv = ures[2,:,:]

# compute udot and vdot (for contour lines). 
# Requires FitzHugh_Nagumo_ODE to be in scope, from 11.6_produce.jl
du = copy(ures)
FitzHugh_Nagumo_ODE(du, ures, params, 0)
udot = du[1,:,:]
vdot = du[2,:,:]

phimat, xps, yps = compute_phases(uu,vv,xx,yy)

# %% plot
fig = plt.figure(figsize=(figx,figy))

ucmap = "gist_earth"
wcmap = "inferno"
pcmap = "hsv"
color_udot = COLORS[1]  # "w"   # 4 o.k.
color_vdot = COLORS[2]   #  6 o.k.

# plot u
axu = subplot(1,4,1)
axu.clear()

im1 = pcolormesh(xx, yy, uu; cmap = ucmap) #"inferno")
axu.contour(xx,yy,udot,levels=[0.],colors="0.5",linewidths=[3.], linestyles = "-")   # with cyan
axu.contour(xx,yy,udot,levels=[0.],colors=color_udot,linewidths=[2.], linestyles = "-")   # with cyan
im1.set_clim(-0.2, 1.)

axu.set_xticks([]) 
axu.set_yticks([]) 
axu.set_aspect("equal")

cbar = colorbar(im1; ax = axu, ticks=[-0.2, 1.], orientation="horizontal", pad=0.02)
cbar.set_label(L"u", labelpad = -30) 
cbar.ax.set_xticklabels(["-0.2", "1"])  # horizontal colorbar

# plot v
axv = subplot(1,4,2)
im2 = axv.pcolormesh(xx,yy,vv,cmap = "inferno")
axv.contour(xx,yy,vdot,levels=[0.],colors=color_vdot,linewidths=[1.5])
im2.set_clim(0., 0.3)

axv.set_xticks([]) 
axv.set_yticks([]) 
axv.set_aspect("equal")

cbar = colorbar(im2; ax = axv, ticks=[0, 0.3], orientation="horizontal", pad=0.02)
cbar.set_label(L"w", labelpad = -30) 
cbar.ax.set_xticklabels(["0", "0.3"])  # horizontal colorbar

# phases
axp = subplot(1,4,3)
axz = subplot(1,4,4)

im3 = nothing

for (j, ax) in enumerate((axp, axz))
    ax.clear()
    global im3 = ax.pcolormesh(xx,yy,phimat,cmap = pcmap) #  cyclic color map
    ax.contour(xx,yy,udot,levels=[0.],colors="0.5",linewidths=[3.], linestyles = "-")   # with cyan
    ax.contour(xx,yy,udot,levels=[0.],colors=color_udot,linewidths=[2.], linestyles = "-")   # with cyan
    ax.contour(xx,yy,vdot,levels=[0.],colors=color_vdot,linewidths=[3])
    im3.set_clim(-pi, pi)

    # add markers at location(s) of PS
    # --------------------------------
    for ips = 1:length(yps)
        ax.scatter(yps[ips],xps[ips]; marker = "o" , s=100, 
        c="0.5", linewidths=2,  zorder = 2, edgecolors = "w")
    end
        
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    ax.set_aspect("equal")

    cbar = colorbar(im3; ax=ax, ticks=[-pi, pi], orientation="horizontal", pad=0.02)
    cbar.ax.set_xticklabels([L"-\pi", L"\pi"])  # horizontal colorbar
    cbar.set_label(L"\alpha", labelpad=-30)

    if j == 1
        # add colorbar
    end
end

dxi = 40
xmin = xx[115] # 220 # 
xmax = xx[115+dxi] # 300 # N = 600
ymin = yy[157] # 320 # N = 600
ymax = yy[157+dxi] # 400 # N = 600

zbox = ((xmin, xmax), (ymin, ymax))


axis_zoomin!(axz, axp, zbox, zbox, "k"; lw = 3)


add_identifiers!(fig, (axu, axv, axp, axz))

fig.subplots_adjust(wspace = 0.15, hspace = 0.51, left = 0.025, 
right = 0.975, bottom = 0, top = 1.5)

wsave(plotsdir("11", "FitzHugh_Nagumo_PS"), fig)
