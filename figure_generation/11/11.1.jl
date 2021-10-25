# Fig_1_Brusselator_eigenvalues_v3.jl
# ==================================== u.p.  23.8.21

using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, OrdinaryDiffEq
using StaticArrays
using LinearAlgebra


# parameters
# ----------
a = 9
bs = [9.8, 10.2]
dvecs = [
    [4.5, 3.0, 2.0, 1.7],
    reverse([0.1, 1.0, 1.7, 2.0]),
]
fig, axs = subplots(1,2)

for (j, b) in enumerate(bs)

dvec = dvecs[j]
ivec = [1 2 6 4]   # order of colors

qmin = 0.0
qmax = 2.0
nq = 301
qvec = range(qmin,qmax; length=nq)
Re_lambda_vec = zeros(nq)
Im_lambda_vec = zeros(nq)
Re_qvec = zeros(nq)
Im_qvec = zeros(nq)

ax = axs[j]

ax.axhline(0; linewidth=1.0,color="k" )

for id = 1:length(dvec)

    dd = dvec[id]
    jj = kk = 0
    for iq = 1:nq
        cc = qvec[iq] * qvec[iq]
        b11 = b - 1 - cc
        b12 = a
        b21 = -b
        b22 = - a - dd*cc
        traceB = b11+b22
        detB = b11*b22 - b12*b21
        rootarg = traceB*traceB/4 - detB
        if rootarg < 0
            Re_lambda = traceB/2
            jj = jj + 1
            Im_qvec[jj] = qvec[iq]
            Im_lambda_vec[jj] = Re_lambda
        else
            Re_lambda = traceB/2 + sqrt(rootarg)
            if kk == 0
               kk = 1
               Re_qvec[kk] = Im_qvec[jj]
               Re_lambda_vec[kk] = Im_lambda_vec[jj]
            end
            kk = kk + 1
            Re_qvec[kk] = qvec[iq]
            Re_lambda_vec[kk] = Re_lambda
        end
    end
    ax.plot(Re_qvec[1:kk],Re_lambda_vec[1:kk], linewidth=3.0,linestyle="--",color=COLORS[ivec[id]] )
    ax.plot(Im_qvec[1:jj],Im_lambda_vec[1:jj], linewidth=3.0,linestyle="-",color=COLORS[ivec[id]] )
end


ax.set_xlim([qmin,qmax])
if j == 1
    ax.set_ylim([-2.1,2.6])
else
    ax.set_ylim([-2.1,0.4])
end
ax.legend([ "_" ; "_" ; L"d="*string(dvec[1]);
                  "_" ; L"d="*string(dvec[2]);
                  "_" ; L"d="*string(dvec[3]);
                  "_" ; L"d="*string(dvec[4]);
              #    "_" ; L"d="*string(dvec[5]);
                  ], 
    fontsize = 26, handlelength = 1, handletextpad = 0.6,
    # loc = "upper left", 
)

end

axs[1].set_xlabel(L"q") 
axs[2].set_xlabel(L"q") 
axs[1].set_ylabel(L"\max (\rm{Re} (\lambda))") 

add_identifiers!(fig)
fig.tight_layout(pad=0.3)
wsave(plotsdir("11", "brusselator_eigenvalues"), fig)
