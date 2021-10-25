# Noise radius illustration using e.g. Poincare section
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

ro = Systems.roessler([0.1,0.2,0.1])
fig, axs = subplots(2, 1; figsize = (0.4*figx, 1.5figy))
plane = (2, 0.0)

err = (1e-4, 1e-12)
εs = 10 .^ (-5:0.5:1)

C = ["C0", "C2"]

for (i, e) ∈ enumerate(err)
    p = poincaresos(
        ro, plane, 10000; Ttr = 100.0,
        rootkw = (xrtol = e, atol = e)
    )

    axs[1].scatter(p[:, 1], p[:, 3]; s = 20, alpha = 0.75/i^2,
    label = "tol = \$10^{$(round(Int, log10(e)))}\$", color = C[i])
    Cs = correlationsum(p, εs)
    y = log.(Cs ./ maximum(Cs))
    axs[2].plot(log.(εs), y; color = C[i])

    # add noise line
    if i == 1
        z = log.(εs)[end-3]
        axs[2].plot([z, z], [-14, y[end-3]]; color = "C1", lw = 2, ls = ":")
        axs[2].text(z+0.2, -9, "\$\\sigma\$"; color = "C1")
    end
end

axs[1].set_xlabel("\$x\$";labelpad = -15)
axs[1].set_ylabel("\$z\$")
axs[1].legend(markerscale = 5, fontsize=26)
axs[1].set_yticks(0:6:18)

axs[2].plot([-9, 0], 0.9 .* [-9, 0]; ls = "--", color = C[2])
axs[2].text(-7, -4, "0.9"; color = C[2])
axs[2].plot([-6, -3], 2.0 .* [-6, -3] .- 2; ls = "--", color = C[1])
axs[2].text(-3, -11, "2"; color = C[1])
axs[2].set_ylabel("\$\\log(C)\$"; labelpad = -10)
axs[2].set_xlabel("\$\\log(\\varepsilon)\$"; labelpad = -12)
axs[2].set_xticks(-12:5:2)
axs[2].set_yticks(-14:5:2)

fig.tight_layout(pad = 0.5)
add_identifiers!(fig)
wsave(plotsdir("5", "fractaldim_noise"), fig)
