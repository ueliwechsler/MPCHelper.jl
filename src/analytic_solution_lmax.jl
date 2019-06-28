# Functionalities for analytical Solution for lmax assuming we have a linear cost
# function ℓ, a Ellispoidal (Convec) set Ω and linear Feedback gain K

using LazySets

"""
        get_Ellipsoid(P, α)

    Creates Ellipsoid `Ω := {x |x'Px ≤ α}` using LazySets"""
function get_Ellipsoid(P, α)
    if !issymmetric(P)
        # if the matrix is not Symmetric at all
        (norm(P-P') >= 1e-10) &&  throw(DomainError(:P, "Matrix needs to be symmetric."))
        @warn "Symmetric(P) is called on P, which is symmetric with numerical noise."
        # if there is numerical noise
    end
    # P_shape = αP⁻¹
    P_shape = α*Symmetric(P)^(-1)
    return Ellipsoid(P_shape)
end

"""
        ω_opt(γ, P, α)

    Analytical Solution for support vector `ω* = arg max γ'ω s.t. ω ∈ Ω`
    for an ellipsoidal set `Ω = {x ∈ ℝⁿ | x'Px ≤ α}` """
ω_opt(γ, P, α) = P^(-1)*γ/norm(P^(-1/2)*γ)*sqrt(α)

"""
        hΩ_opt(γ, P, α)

    Analytical Solution for support function `hΩ(γ) = max γ'ω s.t. ω ∈ Ω`
    for an ellipsoidal set `Ω = {x ∈ ℝⁿ | x'Px ≤ α}` """
hΩ_opt(γ, P, α) = norm(P^(-1/2)*γ)*sqrt(α)

"""
        lmax_mask(z, v, P, α, l_x, l_u, K)

    Calculate `lmax(z,v) = lmax(z, v, P, α, l_x, l_u, K)` for an ellipsoidal set
    `Ω = {x ∈ ℝⁿ | x'Px ≤ α}`, the linear system `l(x,u)= lxᵀx + luᵀu` and the
    linear feedback gain `K`.
    `lmax(z,v) = l(z,v) +  max (lx + K'lu)'ω s.t. ω ∈ Ω`. """
function lmax_mask(z, v, P::AbstractMatrix, α::Real, l_x, l_u, K)
    l(x,u) = dot(l_x, x) + dot(l_u, u)
    γ = l_x + K'*l_u
    hΩ = hΩ_opt(γ, P, α)
    return l(z,v) + hΩ
    # ω = ω_opt(γ, P, α)
    # return l(z,v) + l(ω, K*ω)
end

function lmax_mask(z, v, Ω::LazySet, l_x, l_u, K)
    l(x,u) = dot(l_x, x) + dot(l_u, u)
    γ = l_x + K'*l_u
    hΩ = ρ(vec(γ), Ω)
    return l(z,v) + hΩ
end


# untested max_cost_diff_over_set
"""
        max_cost_diff_over_set(P, α, l_x, l_u, K)

    Calculate `c := max lmax(z, v) - l(z + ω,Kω + v) s.t. (z,v) ∈ ℤ_bar, ω ∈ Ω.`
    For an ellipsoidal set Ω = {x ∈ ℝⁿ | x'Px ≤ α},
    a linear system l(x,u)= lx'x + lu'u and a linear feedback gain K.
    There exists an analytical solution:
    `c = hΩ(lx+Kᵀlu) + hΩ(-(lx+Kᵀlu)) = 2|P^(-1/2)(lx+Kᵀlu)|√(α)`. """
function max_cost_diff_over_set(P::AbstractMatrix, α::Real, l_x, l_u, K)
    l(x,u) = dot(l_x, x) + dot(l_u, u)
    γ = l_x + K'*l_u
    return 2*hΩ_opt(γ, P, α)
    # hΩmax = hΩ_opt(γ, P, α)
    # hΩmin = hΩ_opt(-γ, P, α)
    # return hΩmax + hΩmin
end

function max_cost_diff_over_set(Ω::LazySet, l_x, l_u, K)
    l(x,u) = dot(l_x, x) + dot(l_u, u)
    γ = l_x + K'*l_u
    hΩmax = ρ(vec(γ), Ω)
    hΩmin = ρ(vec(-γ), Ω)
    return hΩmax + hΩmin
end
