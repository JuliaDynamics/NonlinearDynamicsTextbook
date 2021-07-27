using LinearAlgebra: norm, inv, eigvals
import ForwardDiff
using StaticArrays: setindex
using DynamicalSystems

# This code is written for `SVector`s, or out-of-place dynamical systems.
# Notation z is used for the vector (x..., p), i.e. mixed space
"""
    bifurcation_rule_form(ds) → f, J
Produce dynamic rule `f` and Jacobian function `J` in a two-argument form.
"""
function bifurcation_rule_form(ds::ContinuousDynamicalSystem)
    f = (x, p) -> ds.f(x, p, 0.0)
    J = (x, p) -> ds.jacobian(x, p, 0.0)
    return f, J
end

function bifurcation_rule_form(ds::DiscreteDynamicalSystem)
    error("TODO: Convert `f` to `f - x`.")
    f = (x, p) -> ds.f(x, p, 0.0)
    J = (x, p) -> ds.jacobian(x, p, 0.0)
    return f, J
end


# Corrector: Newton's method in mixed-space with extra fixed condition
"""
    corrector(zpred, zprev, f, J; max_steps = 200, δ = 0.9, ε = 1e-3) → (z*, success)
Implement the Newton algorithm to converge to the fixed point for given predicted
mixed state `zpred`, previously found state `zprev`, function `f`, Jacobian `J` 
and stepping factor `δ`.

If iteration occurs for more than `max_steps` without convergence better than `ε`,
then iteration stops and `success` will be `false`.
"""
function corrector(zpred, zprev, f, J; δ = 0.9, max_steps = 200, ε = 1e-6)
    Δz = zpred .- zprev
    # index of variable that changed the most
    # The `iszero` clause is just to keep the first iteration at given `p0`
    i = all(iszero, Δz) ? length(Δz) : argmax(Δz .^ 2)
    c = 0
    zⱼ = zpred
    zⱼ₊₁ = newton_step!(zⱼ, zpred, i, f, J, δ)
    while norm(zⱼ₊₁ - zⱼ) > ε
        zⱼ = zⱼ₊₁
        zⱼ₊₁ = newton_step!(zⱼ, zpred, i, f, J, δ)
        c += 1
        if c > max_steps 
            @warn("Newton did not converge.")
            return (zⱼ₊₁, false)
        end
    end
    return zⱼ₊₁, true
end

function newton_step!(zⱼ, zpred, i, f, J, δ)
    Jfinal = mixed_jacobian(zⱼ, i, f, J)
    xⱼ = zⱼ[1:end-1]; pⱼ = zⱼ[end]
    g = f(xⱼ, pⱼ)
    gz = vcat(g, zⱼ[i] - zpred[i])
    zⱼ₊₁ = zⱼ - δ*(inv(Jfinal))*gz
    return zⱼ₊₁
end

function mixed_jacobian(z, i, f, J)
    x = z[1:end-1]; p = z[end]
    # start creating the mixed space jacobian
    j = J(x, p)
    # to the state space jacobian add one more column, derivative towards p
    pder = ForwardDiff.derivative(p -> f(x, p), p)
    Jmixed = hcat(j, pder)
    # add the last row, which is 1 for the `i` entry, 0 everywhere else
    last_row = setindex((@SVector zeros(length(z))), 1.0, i)
    Jfinal = vcat(Jmixed, last_row')
    return Jfinal
end

# Predictor: Secant (linear extrapolation of already found bifurcation curve)
function predictor(zs, dz0)
    if length(zs) == 1
        return zs[end]
    elseif length(zs) == 2 # 1 entry is z0, 2nd entry is 1st found fixed point
        return zs[end] .+ dz0
    else
        return 2zs[end] .- zs[end-1]
    end
end

# Continuation function: perform a step of predictor-corrector and save values
function continuation!(zs, f, J; dz0, pmin, pmax)
    zpred = predictor(zs, dz0)
    (pmin ≤ zpred[end] ≤ pmax) || return false
    zˣ, success = corrector(zpred, zs[end], f, J)
    push!(zs, zˣ)
    return success
end

# Continuation loop: do continuation for a given amount of steps
function continuation(f, J, x0, p0;
        pmin, pmax, dp0, dx0,
    )

    z0 = vcat(x0, p0); zs = [z0]; dz0 = vcat(dx0, dp0)

    ps = [p0]
    xs = Dataset([x0])
    stability = Bool[]
    for i in 1:N
        success = continuation!(zs, f, J; dz0, pmin, pmax)
        # Stop iteration if we exceed given parameter margins
        success || break
        # Detect stability of found fixed point (needs `Array` coz of StaticArrays.jl)
        eigenvalues = eigvals(Array(J(zs[end][1:end-1], zs[end][end])))
        μ = maximum(real(v) for v ∈ eigenvalues)
        isstable = μ < 0
        push!(stability, isstable)
    end
    xs = Dataset([z[1:end-1] for z in zs])
    ps = [z[end] for z in zs]
    popfirst!(xs.data); popfirst!(ps) # remove initial guess
    return xs, ps, stability
end

