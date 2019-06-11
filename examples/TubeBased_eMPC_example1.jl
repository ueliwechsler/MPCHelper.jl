# Example 1 of Tube Based EMPC Paper
#TODO clean up and add compute_tightened_sets to library!

using MPCHelper

using JuMP
using Plots
using Ipopt
using LinearAlgebra
using LazySets

## Initialize System
A = [0.7776 -0.0045; 26.6185 1.8555]
B = [-0.0004 0.2907]'
W = [-0.0002 0.0893; 0.139 1.2267]
K = [-0.6457 -5.4157]

Aʷ = [1 0.0; 0 1; -1 0; 0 -1]
bʷ = [0.1; 0.1; 0.1; 0.1]
 = linear_map(W, HPolyhedron(Aʷ, bʷ)) # corresponds to W⁻¹Aʷ <= bʷ
plot(𝒲, 1e-3, alpha=0.1)
plot!(𝒲, 1e-3, alpha=0.1)
a = 𝒲 + 𝒲
plot(𝒲, 1e-5)


Aˣ = [1 0; -1 0 ; 0 1.0; 0 -1.0]
bˣ = [0.5; 0.5; 5; 5]
𝕏 = HPolyhedron(Aˣ, bˣ)

# Aᵘ = [1 0; -1 0]
# bᵘ = [15, 15]
# 𝕌 = HPolyhedron(Aᵘ, bᵘ)
𝕌 = Interval(-15.0, 15.0)

# Set tightening for Linear System 2D
function compute_tightened_sets(A, B, K, 𝒲, 𝕏, 𝕌::HPolyhedron; err_approx=5e-6, ϵ=1e-5)
    A_stab = A + B*K  #eigvals(A_stab)
    Ω = approx_mRPI(ϵ, A_stab, 𝒲, err_approx=err_approx)
    𝕏_bar = 𝕏 - Ω
    KΩ = overapproximate(K*Ω, err_approx)
    𝕌_bar = 𝕌 - KΩ
    return 𝕏_bar, 𝕌_bar,  Ω
end

# if 𝕌 is only 1D, i.e. a interval!
function compute_tightened_sets(A, B, K, 𝒲, 𝕏, 𝕌; err_approx=5e-6, ϵ=1e-5)
    A_stab = A + B*K  #eigvals(A_stab)
    Ω = approx_mRPI(ϵ, A_stab, 𝒲, err_approx=err_approx)
    𝕏_bar = 𝕏 - Ω
    KΩ = linear_map(K,Ω)
    𝕌_bar = 𝕌 - KΩ
    !(𝕌_bar ⊆ 𝕌) && error("Interval calculation of 𝕌_bar went wrong!")
    return 𝕏_bar, 𝕌_bar,  Ω
end


@time 𝕏_bar, 𝕌_bar, Ω = compute_tightened_sets(A, B, K, 𝒲, 𝕏, 𝕌)

function l(x1,x2)
    l = x1
    if x2 < -2
        l += 10*x2^2 + 40*x2 + 40
    elseif -2 <= x2 <= 2
        l += 0
    else
        l += 10*x2^2 - 40*x2 + 40
    end
end

## Setup Model

N = 20
n, m = size(B)
z_0 = [-0.258,  3.5]
z_0 = [0,  0]

z_sr = [-0.0076, 0.4275]
z_s = [-0.0358, 2.0009]
z_s = [-0.0675689, 3.78769]
##
mpc_run = Model(with_optimizer(Ipopt.Optimizer))
@variable(mpc_run, z[1:n, 0:N])
@variable(mpc_run, u[1:m, 0:N-1])
constrain_variable_by_set!(mpc_run, z, 𝕏_bar)
constrain_variable_by_set!(mpc_run, u, 𝕌_bar)

# Add dynamics constraints
@constraint(mpc_run, [i=0:N-1, j=1:n],
            z[j, i+1] == A[j,1]*z[1, i] + A[j,2]*z[2, i] + B[j]*u[1,i] )

end_con =     @constraint(mpc_run, [i=1:n], z[i,N] == z_s[i])
initial_con = @constraint(mpc_run, [i=1:n], z[i,0] == z_0[i])

JuMP.register(mpc_run, :l, 2, l, autodiff=true)
@NLobjective(mpc_run, Min, sum(l(z[1,i] - z_s[1], z[2,i] - z_s[2]) for i=0:N) )
# @NLobjective(mpc_run, Min, sum(l(z[1,i], z[2,i]) for i=0:N) )
optimize!(mpc_run)
objective_value(mpc_run)
##

z_real = value.(z).data
u_real = value.(u).data

## Visualize
using Plots
norm_1 = maximum(abs.(z_real[1,:]))
norm_2 = maximum(abs.(z_real[2,:]))
plot(z_real[1,:] , z_real[2,:])
plot!(u_real[1,:])


anim = @gif for i = 1:N
    plot(𝕏_bar, 1e-3, alpha=0.1)
    plot!(z_real[1, :], z_real[2, :], xlim=(-0.5, 0.5), ylim=(-5, 5))
    plot!([z_real[1, i]], [z_real[2, i]], marker=(:hex, 6))
end


list_of_constraint_types(mpc_run)
all_constraints(opt, GenericAffExpr{Float64, VariableRef}, MOI.EqualTo{Float64})
a = all_constraints(mpc_run, GenericAffExpr{Float64, VariableRef}, MOI.LessThan{Float64})





## Animate MPC
z_s = [-0.0675689, 3.78769]
function run_mpc(z_0)
    mpc_run = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(mpc_run, z[1:n, 0:N])
    @variable(mpc_run, v[1:m, 0:N-1])
    constrain_variable_by_set!(mpc_run, z, 𝕏_bar)
    constrain_variable_by_set!(mpc_run, v, 𝕌_bar)

    # Add dynamics constraints
    @constraint(mpc_run, [i=0:N-1, j=1:n],
                z[j, i+1] == A[j,1]*z[1, i] + A[j,2]*z[2, i] + B[j]*v[1,i] )

    end_con =     @constraint(mpc_run, [i=1:n], z[i,N] == z_s[i])
    initial_con = @constraint(mpc_run, [i=1:n], z[i,0] == z_0[i])

    JuMP.register(mpc_run, :l, 2, l, autodiff=true)
    @NLobjective(mpc_run, Min, sum(l(z[1,i] - z_s[1], z[2,i] - z_s[2]) for i=0:N) )
    # @NLobjective(mpc_run, Min, sum(l(z[1,i], z[2,i]) for i=0:N) )
    optimize!(mpc_run)
    # objective_value(mpc_run)
    return value.(z).data, value.(v).data
end


# The robot's starting position and velocity
function does()
z= [0,  0]
x = [0,  0]
@gif for i in 1:20
    @show z
    # Plot the current position
    plot(z[1, :], z[2, :], marker=(:hex, 4), xlim=(-0.5, 0.5), ylim=(-5, 5))
    # Run the MPC control optimization
    z_plan, v = run_mpc(z)

    # Draw the planned future states from the MPC optimization
    plot!(z_plan[1, :], z_plan[2, :], linewidth=2)

    # Apply the planned acceleration and simulate one step in time
    v = v[:, 1]
    z⁺ = A*z + B*v
    u = K*(x-z) + v
    x⁺ = A*x + B*u
    @show z⁺, x⁺
    z = z⁺
    x = x⁺
end
end

function simulate_system(A, B, K, 𝒲)


end

does()


plot(𝕏 - Ω + Ω, 1e-3)
plot!(𝕏 - Ω, 1e-3)
