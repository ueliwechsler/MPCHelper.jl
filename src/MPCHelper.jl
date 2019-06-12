module MPCHelper

using JuMP
using LinearAlgebra
using LazySets
using Polyhedra
import LazySets.Approximations.overapproximate
import LazySets: Interval
# using ControlSystems
using Plots

include("JuMPHelper.jl")
include("approximate_mRPI_set.jl")
include("set_calculations_for_MPC.jl")
include("analytic_solution_lmax.jl")
include("PlotHelper.jl")

# copied functionalities
# include("ControlSystem.jl")


# approximate_mRPI_set
export approx_mRPI
# tightend_sets
export compute_tightened_sets
# JuMPHelper
export constrain_variable_by_set!
# PlotHelper
export plot_lin_contour!, plot_vec!
# analytic_solution_lmax
export lmax_mask, ω_opt, hΩ_opt, get_Ellipsoid


end # module
