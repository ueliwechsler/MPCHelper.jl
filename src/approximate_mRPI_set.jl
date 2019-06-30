## Implementation of the algorithm described in: "Invariant approximation of the minimal robust positively invariant set" S.V.Rakovic
# Invariant approximations of the minimal robust positively invariant set for and2D systems
#TODO 1D-systems

using LinearAlgebra
using LazySets
using Polyhedra
using LazySets.Approximations # overapproximate, SphericalDirections

""" Pontryagin Set Difference for HPolytope and HPolygon
    From Theory and Computation of Disturbance Invariant Sets for discrete time
    linear systems https://www.hindawi.com/journals/mpe/1998/934097/abs/"""
Base.:-(X::HPolytope, Y::HPolytope) = begin
    ax, bx, n = get_constraints(X)
    b_xminusy = bx .- ρ.(ax,Ref(Y))
    XminusY = HPolytope(vcat(ax'...), b_xminusy)
    return XminusY
end

Base.:-(X::AbstractPolytope, Y::AbstractPolytope) = begin
    ax, bx, n = get_constraints(X)
    b_xminusy = bx .- ρ.(ax,Ref(Y))
    XminusY = HPolytope(vcat(ax'...), b_xminusy)
    return XminusY
end


""" Pontryagin Set Difference for two Interval from IntervalArithmetic"""
Base.:-(a::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}},
        b::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}}) = begin
    low = a.dat.lo - b.dat.lo
    high = a.dat.hi - b.dat.hi
    return LazySets.Interval(low, high)
end

function get_support(P::AbstractPolytope, dims)
    x = Vector{Tuple{Float64,Float64}}(undef,dims)
    E = Matrix(I,dims,dims)
    for i=1:dims
        x[i] = (-ρ(-e(i,dims), P), ρ(e(i,dims), P))
    end
    return x
end

function get_constraints(𝒫::AbstractPolyhedron)
    # not all AbstractPolyhedron have 𝒫.constraints => constraints_list
    a = getfield.(constraints_list(𝒫),:a) # normal vector of halfspace
    b = getfield.(constraints_list(𝒫),:b) # support value of halfspace
    n = length(𝒫.constraints)
    return a,b, n
end

""" Equation (11) """
function αᴼ(s::Int64, A::AbstractMatrix, 𝒫::HPolytope)
    a, b, n_halfspaces = get_constraints(𝒫)
    hw_vec = zeros(n_halfspaces)
    for i=1:n_halfspaces
        hw_vec[i] = ρ((A^s)'*a[i], 𝒫) / b[i]
    end
    return maximum(hw_vec)
end

## Approximation to Invariant Set
"""  Equation (13) """
function M(s::Int64, A::AbstractMatrix, 𝒫::HPolytope)
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

# linear_map(A, P) is efficient, if A is invertible, since it does not have to convert
# to V-Representation!
function F(s::Int64, A::AbstractMatrix, 𝒫::HPolytope)
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

#TODO α and M for Interval type (1D)
# SinglePoint: A = Interval(0.5,0.5)
# or Apply just geometric series 1/(1-\rho)
function F(s::Int64, A::Real, 𝒫::Interval)
    Aint = Interval(A,A)
    # As defined in the paper, but different loop
    if s == 0
        return Interval(0.0,0.0)
    else
        F = 𝒫
        Ai = Aint
        for i=1:s-1
            F += 𝒫*Ai
            Ai = Ai*Aint
        end
        return F
    end
end

# The overapproximate case only works in the case of a HPolygon!
"""
    approx_mRPI(ϵ::Float64, A::AbstractMatrix, 𝒲::HPolyhedron;
                err_approx::Float64=1e-5, isLazy::Bool=false)

Calculate approximative mRPI set according to "Invariant Approximation of the mRPI set" by Rakovic et all."""
function approx_mRPI(ϵ::Float64, A::AbstractMatrix, 𝒲::HPolygon;
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
    return overapproximate(lazyF_infty, err_approx) #overapproximate only works for 2D
end

# Only works for 3D 🙈
function approx_mRPI(ϵ::Float64, A::AbstractMatrix, 𝒲::HPolytope;
                    n_dirs::Int=10, isLazy::Bool=false)
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
    dirs = SphericalDirections(n_dirs, n_dirs)
    return overapproximate(lazyF_infty, dirs)
end


function approx_mRPI(s::Int, A::Real, 𝒲::Interval)
    return F(s, A, 𝒲)
end

## use https://github.com/JuliaPolyhedra/Polyhedra.jl/blob/master/examples/Minimal%20Robust%20Positively%20Invariant%20Set.ipynb for 3D and more!

# W = polyhedron(Wv, CDDLib.Library())

# using Polyhedra
# Wv = vrep([[x, y] for x in [-1.0, 1.0] for y in [-1.0, 1.0]])
# W = polyhedron(Wv)
#
# function vertices_to_polhedron(x_vec,y_vec)
#     Wv = vrep([[x, y] for x in x_vec for y in y_vec])
#     W = polyhedron(Wv)
# end
#
# function vertices_to_polhedron(x_vec,y_vec, z_vec)
#     Wv = vrep([[x, y,z] for x in x_vec for y in y_vec for z in z_vec])
#     W = polyhedron(Wv)
# end
#
# function Fs(s::Integer, verbose=1)
#     @assert s ≥ 1
#     F = W
#     A_W = W
#     for i in 1:(s-1)
#         A_W = A * A_W
#         F += A_W
#         if verbose ≥ 1
#             println("Number of points after adding A^$i * W: ", npoints(F))
#         end
#         removevredundancy!(F)
#         if verbose ≥ 1
#             println("Number of points after removing redundant ones: ", npoints(F))
#         end
#     end
#     return F
# end
