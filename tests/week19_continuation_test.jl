include("../examples/example_functions.jl")
include("../numerical_continuation.jl")


# Excercise 1
algebraic_pv_1, algebraic_c_1 = continuation(algebraic, [-1 0], 6, "c", -2:0.01:2,
    method="natural_parameter_continuation", discretisation="none")

hopf2d_pv_1, hopf2d_c_1 = continuation(hopf2d, [-1 0], 6, "beta", 0:0.01:2,
    method="natural_parameter_continuation", discretisation="shooting")

hopf2d_mod_pv_1, hopf2d_mod_c_1 = continuation(hopf2d_modified, [-1 0], 6, "beta", 0:0.1:2,
    method="natural_parameter_continuation", discretisation="shooting")


# Excercise 2
algebraic_pv_2, algebraic_c_2 = continuation(algebraic, [1 1], 6, "c", -2:0.01:2,
    method="pseudo_arclength", discretisation="none")

hopf2d_pv_2, hopf2d_c_2 = continuation(hopf2d, [-1 0], 6, "beta", 0:0.01:2,
    method="pseudo_arclength", discretisation="shooting")

hopf2d_mod_pv_2, hopf2d_mod_c_2 = continuation(hopf2d_modified, [-1 0], 6, "beta", 0:0.1:2,
    method="pseudo_arclength", discretisation="shooting")
