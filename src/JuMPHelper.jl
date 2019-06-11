# using JuMP


# for n>1 nD constraint set 𝒳 and nD Variable x
function constrain_variable_by_set!(model::Model, x::AbstractArray{VariableRef,2}, 𝒳::AbstractPolyhedron)
    a, b, s = get_constraints(𝒳)
    N = size(x,2) - 1
    for j=0:N
        @constraint(model, [i=1:s], a[i]'*x[:, j] <= b[i])
    end
end

function constrain_variable_by_set!(model::Model, x::AbstractArray{VariableRef},
                                    𝒳::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}})
    N = size(x,2)
    a = [1.0, -1.0]
    b = [𝒳.dat.hi, -𝒳.dat.lo]
    for j=0:N-1
        @constraint(model, [i=1:2], a[i]'*x[1, j] <= b[i])
    end
end

# 1D variable x and 1D Interval set 𝒳
function constrain_variable_by_set!(model::Model, x::VariableRef,
                                    𝒳::Interval{Float64,LazySets.IntervalArithmetic.Interval{Float64}})
    a = [1.0, -1.0]
    b = [𝒳.dat.hi, -𝒳.dat.lo]
    @constraint(model, [i=1:2], a[i]'*x <= b[i])
end


function nonlin_opt_okay(model::Model)
    term_status = termination_status(model)
    if !(Int(term_status) == 1 || Int(term_status) == 4)
        error("Optimization not a local or global Optimum: $term_status")
    end
    nothing
end
