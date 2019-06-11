module MPCHelper

using JuMP
using LinearAlgebra
using LazySets
using Polyhedra
import LazySets.Approximations.overapproximate
import LazySets: Interval

include("JuMPHelper.jl")
include("approximate_mRPI_set.jl")
include("tightend_sets.jl")

export approx_mRPI
export compute_tightened_sets
export constrain_variable_by_set!


end # module
