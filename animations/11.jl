using DrWatson
@quickactivate "NonlinearDynamicsTextbook"
include(srcdir("style.jl"))

using DynamicalSystems, PyPlot

# The solution explicitly assumes δL = 1
∂²(u, i, L) = u[mod1(i-1, L)] -2u[i] + u[mod1(i+1, L)]
function ∂⁴(u, i, L)
    u[mod1(i-2, L)] +
    -4*u[mod1(i-1, L)] +
    6u[i] +
    -4u[mod1(i+1, L)] +
    u[mod1(i+2, L)]
end

function swift_hohenberg(T, r, dt, u0)
    L = length(u0)
    U = zeros(T, length(u0))
    U[1, :] .= u0
    for j in 1:T-1
        @show j
        for i in 1:L
            U[j+1, i] = U[j, i] + dt*(
                        (r-1)*U[j, i] -
                        2∂²(view(U, j, :), i, L) -
                        ∂⁴(view(U, j, :), i, L) - U[j, i]^3
            )
        end
    end
    return U
end

dt = 0.05
u0 = rand(100)
r = 0.5
T = 1000

U = swift_hohenberg(T, r, dt, u0)

figure()
imshow(U')
xlabel("time")
ylabel("space")
gca().set_aspect("auto")
colorbar()
tight_layout(;pad = 0.25)

# figure()
# plot(U[end, :])
# xlabel("space")
# ylabel("u")
# tight_layout(;pad = 0.25)
