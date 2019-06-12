# Functionalities for analytical Solution for lmax assuming we have a linear cost
# function ℓ, a Ellispoidal (Convec) set Ω and linear Feedback gain K

using LazySets

# lmax(z,v) = lmax(z,v, P, α, l_x, l_u, K)
# lmax(z,v)
"""
     lmax(z, v, P, α, l_x, l_u, K)

     Calculate lmax(z,v) = lmax(z, v, P, α, l_x, l_u, K) for an ellipsoidal set
     Ω = {x ∈ ℝⁿ | x'Px ≤ α}, the linear system l(x,u)= lx'x + lu'u and the
     linear feedback law K.
     lmax(z,v) = l(z,v) +  max (lx + K'lu)'ω s.t. ω ∈ Ω
     """
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



"""
    ω_opt(γ, P, α)

    Analytical Solution for support vector ω* = arg max γ'ω s.t. ω ∈ Ω
    for an ellipsoidal set Ω = {x ∈ ℝⁿ | x'Px ≤ α} """
ω_opt(γ, P, α) = P^(-1)*γ/norm(P^(-1/2)*γ)*sqrt(α)

"""
    hΩ_opt(γ, P, α)

    Analytical Solution for support function hΩ(γ) = max γ'ω s.t. ω ∈ Ω
    for an ellipsoidal set Ω = {x ∈ ℝⁿ | x'Px ≤ α} """
hΩ_opt(γ, P, α) = norm(P^(-1/2)*γ)*sqrt(α)

"""
    get_Ellipsoid(P, α)

    Creates Ellipsoid Ω := {x |x'Px ≤ α} using LazySets"""
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
