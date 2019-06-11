module MPCHelper

using JuMP
using LinearAlgebra
using LazySets
using Polyhedra
import LazySets.Approximations.overapproximate
import LazySets: Interval

include("JuMPHelper.jl")
include("approximate_mRPI_set.jl")

export approx_mRPI
export constrain_variable_by_set!


end # module
