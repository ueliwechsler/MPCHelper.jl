# HOME_DIR = dirname(@__DIR__)
# PACKAGE_DIR = dirname(HOME_DIR)
# pushfirst!(LOAD_PATH, PACKAGE_DIR)
using MPCHelper
using LazySets

a = HPolytope([eye(2); -eye(2)], [2,2,2,0.0])
b = HPolytope([eye(2); -eye(2)], [0.5,0.5,0.5,0.3])
c = a-b
get_support(a,2)
get_support(c,2)
