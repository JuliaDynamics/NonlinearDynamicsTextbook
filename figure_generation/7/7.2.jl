# %% transfer entropy
using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot

using TransferEntropy, Random

function ulam(dx, x, p, t)
    f(x) = 2 - x^2;
    ε = p[1];
    N = length(x)
    for i in 1:N
        dx[i] = f(ε*x[mod1(i-1, N)] + (1-ε)*x[i])
    end
end
ds = DiscreteDynamicalSystem(ulam, rand(100), [0.04])

genmeth(r) = VisitationFrequency(RectangularBinning(r))
methods = genmeth.((0.01, 0.1, 0.4))
εs = 0.0:0.01:1.0
tes = [zeros(length(εs), 2) for j in 1:length(methods)]

for (i, ε) in enumerate(εs), (j, meth) in enumerate(methods)
    set_parameter!(ds, 1, ε)
    A = trajectory(ds, 10000; Ttr = 10000)
    x1 = A[:, 1]; x2 = A[:, 2]
    tes[j][i, 1] = transferentropy(x1, x2, meth)
    tes[j][i, 2] = transferentropy(x2, x1, meth)
end


fig = figure(figsize = (6figx/10, figx/3))
rs = (0.01, 0.1, 0.4)

for j in 1:length(methods)
    plot(εs, tes[j][:, 1], color = "C$(j-1)", label = "\$r=$(rs[j])\$")
    plot(εs, tes[j][:, 2], color = "C$(j-1)", ls = "dashed")
end
xlabel("coupling strength \$\\epsilon\$")
ylabel("transfer entropy")
legend(loc = "center")
tight_layout(pad=0.3)
wsave(plotsdir("transfer"), fig)