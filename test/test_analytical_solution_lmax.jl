HOME_DIR = dirname(@__DIR__)
PACKAGE_DIR = dirname(HOME_DIR)
pushfirst!(LOAD_PATH, PACKAGE_DIR)
using MPCHelper


using ControlSystems
using LinearAlgebra
using LazySets

using Test

using Plots

# Initializatoin ---------------------------------------------------------------

## Linear 2D System
A = [1 -2; 0 1]
B= Matrix([1 1]')
R = 0.1*I
Q = Matrix(I,size(A))/10
K = dlqr(A,B,Q,R)
A_stab = A - B*K; eigvals(A_stab)
P = dlyap(A_stab, Q)
l_x = [3,3]
l_u = 1
z = [0.0 ;0.0]; v = 0.0

# # Linear 3D System
# A = [0.8 -2 3; 0 0.5 0; 0 0 0.3]
# B= [1 0; 0 1; 0 0.4]
# R = 0.1*Matrix(I,size(B,2), size(B,2))
# Q = 0.1*Matrix(I,size(A))
# K = dlqr(A,B,Q,R)
# A_stab = A - B*K; eigvals(A_stab)
# P = dlyap(A_stab, Q)
# l_x = [3,3,3]
# l_u = [2,2]
# z = [0.0 ;0.0; 0.0]; v = [0.0, 0.0]

# Main Part
α = 2
E = get_Ellipsoid(P, α)
l(x,u) = dot(l_x,x) + dot(l_u,u)
γ = l_x + K'*l_u
ωstar = ω_opt(γ, P, α)
hΩ = hΩ_opt(γ, P, α)
# Scale γ -> vec_γ, such that vec_γ'*γ the same cost as ω*γ
vec_γ = γ*hΩ/norm(γ)^2
@test ωstar'*γ ≈ vec_γ'*γ
# lmax(z,v) for Ellispoid defined by P, α
lmax(z,v) = lmax_mask(z,v, P, α, l_x, l_u, K)
# lmax(z,v) for Ellispoid defined by Ellipsoid E (Type LazySet)
lmax_lazy(z,v) = lmax_mask(z,v, E, l_x, l_u, K)
lmax_opt = lmax(z,v)
lmax2_opt = lmax_lazy(z,v)
@test hΩ ≈ lmax_opt
x_opt = z .+ ωstar
u_opt = vec(v .+ K*ωstar)
@test l(x_opt, u_opt) ≈ lmax_opt

## Plotting (only for 2D systems)
grid_res = -5:0.2:5
plot(E,1e-3, alpha=0.5, lab="Omega", size=(500,500))
plot_vec!(ωstar, line=(3,:blue), lab="w_opt")
plot_vec!(vec_γ, line=(2,:black, :dash), lab="lx + K'lu")
# plot_vec!([1,1], line=(2,:black, :dash), lab="lx + K'lu")
plot_lin_contour!(l_x, grid_res; color=:bluesreds, line=(1,:dash), colorbar=false)
plot_lin_contour!(γ, grid_res; color=:bluesreds, line=(1))#, contour_labels = true)
l_opt(x) = (-γ[1]*x + hΩ) / γ[2] # γ'x = hα
plot!(grid_res, l_opt.(grid_res), line=(:red), lab="lopt")
xlims!(-5,5)
ylims!(-5,5)



## PLotting Surfaces
grid_res = 0:0.1:10
surface(grid_res,grid_res, (x,y) -> lmax([x,y],0) )
surface!(grid_res,grid_res, (x,y) -> l([x,y],0))
