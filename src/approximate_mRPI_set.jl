## Implementation of the algorithm described in: "Invariant approximation of the minimal robust positively invariant set" S.V.Rakovic
# Invariant approximations of the minimal robust positively invariant set for 1D? and2D systems

# using LinearAlgebra
# using LazySets
# using Polyhedra
# import LazySets.Approximations.overapproximate

""" Pontryagin Set Difference for HPolyhedron and HPolygon
    From Theory and Computation of Disturbance Invariant Sets for discrete time
    linear systems https://www.hindawi.com/journals/mpe/1998/934097/abs/"""
Base.:-(X::HPolyhedron, Y::HPolygon) = begin
    ax, bx, n = get_constraints(X)
    b_xminusy = bx .- ρ.(ax,Ref(Y))
    XminusY = HPolyhedron(vcat(ax'...), b_xminusy)
    return XminusY
end

""" Pontryagin Set Difference for two Interval from IntervalArithmetic"""
Base.:-(a::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}},
        b::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}}) = begin
    low = a.dat.lo - b.dat.lo
    high = a.dat.hi - b.dat.hi
    return LazySets.Interval(low, high)
end

function get_constraints(𝒫::AbstractPolyhedron)
    a = getfield.(𝒫.constraints,:a) # normal vector of halfspace
    b = getfield.(𝒫.constraints,:b) # support value of halfspace
    n = length(𝒫.constraints)
    return a,b, n
end

""" Equation (11) """
function αᴼ(s::Int64, A::AbstractMatrix, 𝒫::HPolyhedron)
    a, b, n_halfspaces = get_constraints(𝒫)
    hw_vec = zeros(n_halfspaces)
    for i=1:n_halfspaces
        hw_vec[i] = ρ((A^s)'*a[i], 𝒫) / b[i]
    end
    return maximum(hw_vec)
end

## Approximation to Invariant Set
"""  Equation (13) """
function M(s::Int64, A::AbstractMatrix, 𝒫::HPolyhedron)
    n = size(A, 2)
    E = Matrix{Float64}(I, n, n)
    hw_vec = zeros(n,2)
    for j=1:n
        eⱼ = E[:,j]
        sup_value_1 = 0
        sup_value_2 = 0
        for i=0:s-1
            sup_value_1 += ρ((A^i)'*eⱼ, 𝒫)
            sup_value_2 += ρ(-(A^i)'*eⱼ, 𝒫)
        end
        hw_vec[j,1] = sup_value_1
        hw_vec[j,2] = sup_value_2
    end
    return maximum(hw_vec)
end

function F(s::Int64, A::AbstractMatrix, 𝒫::HPolyhedron)
    # As defined in the paper, but different loop
    if s == 0
        return Ball2(zeros(size(A,1)), 0.0)
    else
        F = 𝒫
        for i=1:s-1
            F += linear_map(A^i,𝒫)
        end
        return F
    end
end


"""
    approx_mRPI(ϵ::Float64, A::AbstractMatrix, 𝒲::HPolyhedron;
                err_approx::Float64=1e-5, isLazy::Bool=false)

Calculate approximative mRPI set according to "Invariant Approximation of the mRPI set" by Rakovic et all."""
function approx_mRPI(ϵ::Float64, A::AbstractMatrix, 𝒲::HPolyhedron;
                    err_approx::Float64=1e-5, isLazy::Bool=false)
    s = 0; α = 1; m = 1
    while α > ϵ/(ϵ + m)
        s += 1
        α = αᴼ(s, A, 𝒲)
        m = M(s, A, 𝒲)
    end
    Fs = F(s, A, 𝒲)
    @show s, α
    lazyF_infty = (1 - α)*Fs # lazy (symoblic) approximative mRPI-set
    # Return the lazy (symbolic) set or an HPolygon
    (isLazy) && return lazyF_infty
    return overapproximate(lazyF_infty, err_approx)
end
