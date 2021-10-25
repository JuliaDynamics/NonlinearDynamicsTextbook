# kuramoto
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq, Random, Statistics

function plot_circle!(ax, phi_vec,omega_vec, add_colorbar=true; s = 1)
    cmap = "viridis"
    # plot circle
    alpha_vec = range(0, 2. * pi; length = 200)
    x_vec = cos.(alpha_vec)
    y_vec = sin.(alpha_vec)
    ax.plot(x_vec,y_vec, linewidth=1.0,color="k") # unit circle  # linestyle="--",
    
    # plot markers
    u_vec = cos.(phi_vec)
    v_vec = sin.(phi_vec)
    implt = ax.scatter(u_vec, v_vec, marker = "o" , s=70, c=omega_vec, 
    cmap = cmap, zorder = 3, edgecolors = "k", linewidths = 1, vmin = -2, vmax = +2)

    # plot arrows
    D = length(u_vec)
    ax.quiver(zeros(D), zeros(D), u_vec, v_vec, omega_vec;
        # color = COLORS[5],
        cmap = cmap, # color = omega_vec,
        scale = 1.5, scale_units = "xy"
    )

    # Plot mean field   
    mf = sum(exp(im*φ) for φ in phi_vec)/length(phi_vec)
    ax.quiver([0.0], [0.0], [s*real(mf)], [s*imag(mf)];
        color = COLORS[4], width = 0.02,
        scale = 1.5, scale_units = "xy"
    )

    ax.set_xticks([]) 
    ax.set_yticks([]) 
    ax.set_aspect("equal")
    if add_colorbar
        cbar = plt.colorbar(implt; ax, orientation = "horizontal", fraction=0.046, pad=0.04)
        cbar.set_ticks([-2,2])
        cbar.set_label(L"\omega_n", labelpad = -30)
    end

end
    
nosc1 = 50     # number of oscillators for (a), (b), (c)
nosc2 = 2000    # number of oscillators for second curve in (c)
sigma = 1.   # 0.1
T = 100.0
Ttr = 100.0
Δt = 0.1
t = 0:Δt:T
nts = length(t)

# open plot
fig = plt.figure(figsize=(figx,figy))

# plot identifiers first, otherwise they appear with colorbars
ax1 = subplot(1,3,1)
ax2 = subplot(1,3,2)
ax3 = subplot(1,3,3)

Kmin = 0.0
Kmax = 8.0   # 0.8
nK = 81
K_vec = LinRange(Kmin,Kmax,nK)
Rmean_vec = zeros(nK)
Rstd_vec = zeros(nK)

for (nj, nosc) in enumerate((nosc1, nosc2))

    rng = MersenneTwister(1234) 
    omega_vec = sigma .* randn(rng,nosc)
    u0 = zeros(nosc)
    ds = Systems.kuramoto(nosc, u0; K = 0.0, ω = omega_vec)

    for iK = 1:nK
        @show iK
        K = K_vec[iK]
        set_parameter!(ds, :K, K)
        tr = trajectory(ds, T; Ttr, Δt)
        phi_vec = tr[end]

        Rts = zeros(nts)
        for j = 1:nts  # R time series
            Rts[j] = abs( sum( exp.(im .* tr[j]) ) ) ./ nosc
        end
        Rmean_vec[iK] = mean(Rts)  
        Rstd_vec[iK] = std(Rts) 

        if iK == 11 && nj == 1    
            plot_circle!(ax1, phi_vec, omega_vec)
        end
        if iK == 81 && nj == 1  # sync
            plot_circle!(ax2, phi_vec,omega_vec)
        end

    end
    ax3.plot(K_vec,Rmean_vec, linewidth=2.0,color=COLORS[nj]) 
    ax3.fill_between(K_vec, Rmean_vec - Rstd_vec, Rmean_vec + Rstd_vec;
        color = COLORS[nj], alpha = 0.15
    )
end

ax3.set_ylabel("\$R\$")
ax3.set_xticks([0, 2.5, 5, 7.5])
ax3.set_xticklabels(["0", "2.5", "5", "7.5"])
ax3.set_xlabel("\$K\$", labelpad = -20)

add_identifiers!(fig, (ax1,ax2,ax3))  
fig.tight_layout(pad=0.38)
wsave(plotsdir("9", "kuramoto"), fig)
