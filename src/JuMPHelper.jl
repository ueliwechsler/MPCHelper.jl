# using JuMP
# using LazySets
#TODO: You could use Polyhedra to get the feasible set (poly = polyhedron(m) m=Model())
#TODO: Naming of the constraints => in order to have better debugging

# Used for DenseAxisArray ========================================
# x ∈ ℝn×N+1 and 𝒳 ⊂ ℝⁿ
function add_constraint!(m::Model, x::AbstractArray{VariableRef,2}, 𝒳::AbstractPolyhedron)
    a, b, s = get_constraints(𝒳)
    N = size(x,2)-1 # Since the first element of the JumpVariable is 0
    for j=0:N
        @constraint(m, [i=1:s], a[i]'*x[:, j] <= b[i])
    end
end

# x ∈ ℝⁿ and 𝒳 ⊂ ℝⁿ
function add_constraint!(m::Model, x::AbstractVector{VariableRef}, 𝒳::AbstractPolyhedron)
    a, b, s = get_constraints(𝒳)
    @constraint(m, [i=1:s], a[i]'*x <= b[i])
end

# x ∈ ℝⁿ and 𝒳 = 𝓍 ∈ ℝⁿ
function add_constraint!(m::Model, x, 𝓍::AbstractVector)
    n = length(𝓍)
    @constraint(m, [i=1:n], x[i] == 𝓍[i])
end

# Used for Array{VariableRef, 2} ========================================

# x ∈ ℝⁿ and 𝒳 = 𝓍 ∈ ℝⁿ
function add_constraint!(m::Model, x, 𝓍::AbstractVector)
    n = length(𝓍)
    @constraint(m, x .== 𝓍)
end

# # x ∈ ℝⁿ and 𝒳 = 𝓍 ± ϵ ∈ ℝⁿ (relaxed terminal constraint)
# function add_constraint!(m::Model, x, 𝓍::AbstractVector, ϵ=1e-5)
#     n = length(𝓍)
#     @constraint(m, [i=1:n], x[i] <= 𝓍[i] + ϵ)
#     @constraint(m, [i=1:n], x[i] <= 𝓍[i] - ϵ)
# end

# # TODO: make it a macro! such that the constr name can be added to the constraint
# function constrain_variable_by_set!(model::Model, x::AbstractVector{VariableRef},
#                                     𝒳::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}})
#     N = length(x)
#     a = [1.0, -1.0]
#     b = [𝒳.dat.hi, -𝒳.dat.lo]
#     for j=0:N-1
#         @constraint(model, [i=1:2], a[i]'*x[j] <= b[i])
#     end
# end
#
# # 1D variable x and 1D Interval set 𝒳
# function constrain_variable_by_set!(model::Model, x::VariableRef,
#                                     𝒳::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}})
#     a = [1.0, -1.0]
#     b = [𝒳.dat.hi, -𝒳.dat.lo]
#     @constraint(model, [i=1:2], a[i]'*x <= b[i])
# end

""" Checks JuMP.MOI.TerminationStatusCode """
function nonlin_opt_okay(model::Model)
    term_status = termination_status(model)
    if Int(term_status) == 1
        # everything is fine
    elseif Int(term_status) == 4
        @warn("No global Optimum found: $term_status")
    elseif (Int(term_status) == 7 || Int(term_status) == 10)
        @warn("Solved to relaxed tolerances: $term_status")
    else
        error("Optimization not a local or global Optimum: $term_status")
    end
    nothing
end

function get_ineq_constr(model::Model)
    return all_constraints(model, GenericAffExpr{Float64, VariableRef}, MOI.LessThan{Float64})
end

function get_eq_constr(model::Model)
    return all_constraints(model, GenericAffExpr{Float64, VariableRef}, MOI.EqualTo{Float64})
end

# list_of_constraint_types(model)
# TODO: print name of constraints!
# print(model)
# Min f(l(z + -0.200000000000002, v), l(z + 0.200000000000002, v))
# Subject to
#  0.5 z - v == 0.0
#  z <= 4.799999999999998
#  -z <= 4.799999999999998
#  v <= 2.0
#  -v <= 2.0
