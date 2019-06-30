using LinearAlgebra


eye(n) = Matrix{Float64}(I,n,n)
eye(n,m) = Matrix{Float64}(I,n,m)

function e(i, dims)
    E = Matrix{Float64}(I, dims, dims)
    return E[:,i]
end
