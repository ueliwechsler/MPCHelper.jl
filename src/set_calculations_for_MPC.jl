
## Compute tightened Sets using mRPI Ω

# if m=n=w=1, i.e. 1-D Linear system
function compute_tightened_sets(A::Real, B::Real, K::Real, 𝒲::Interval, 𝕏::Interval, 𝕌::Interval; s=100)
    A_stab = A + B*K  #eigvals(A_stab)
    K = Interval(K,K)
    Ω = approx_mRPI(s, A_stab, 𝒲)
    𝕏_bar = 𝕏 - Ω
    KΩ = K*Ω
    𝕌_bar = 𝕌 - KΩ
    !(𝕌_bar ⊆ 𝕌) && error("Interval calculation of 𝕌_bar went wrong!")
    return 𝕏_bar, 𝕌_bar,  Ω
end

# Set tightening for Linear System 2D #TODO check if it is general enough
function compute_tightened_sets(A, B, K, 𝒲, 𝕏, 𝕌::HPolyhedron; err_approx=5e-6, ϵ=1e-5)
    A_stab = A + B*K  #eigvals(A_stab)
    Ω = approx_mRPI(ϵ, A_stab, 𝒲, err_approx=err_approx)
    𝕏_bar = 𝕏 - Ω
    KΩ = overapproximate(K*Ω, err_approx)
    𝕌_bar = 𝕌 - KΩ
    return 𝕏_bar, 𝕌_bar,  Ω
end
