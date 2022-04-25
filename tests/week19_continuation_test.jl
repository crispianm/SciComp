# Import Subsystems
include("../examples/example_functions.jl")
include("../finite_difference.jl")
include("../numerical_continuation.jl")
include("../numerical_shooting.jl")
include("../ode_solver.jl")
include("../visualisation.jl")


# Excercise 1
algebraic_pv_1, algebraic_c_1 = continuation(algebraic,[1 1],6,"c",-2:0.05:2,method = "natural_parameter_continuation",discretisation = "none",);
hopf2d_pv_1, hopf2d_c_1 = continuation(hopf2d,[-1 0],6,"beta",2:-0.05:0,method = "natural_parameter_continuation",discretisation = "shooting",);
hopf2d_mod_pv_1, hopf2d_mod_c_1 = continuation(hopf2d_modified,[-1 0],6,"beta",2:-0.05:-1,method = "natural_parameter_continuation",discretisation = "shooting",);

# Excercise 2
algebraic_pv_2, algebraic_c_2 = continuation(algebraic,[1 1],6,"c",-2:0.05:2,method = "pseudo_arclength",discretisation = "none",);
hopf2d_pv_2, hopf2d_c_2 = continuation(hopf2d,[-1 0],6,"beta",2:-0.05:0,method = "pseudo_arclength",discretisation = "shooting",);
hopf2d_mod_pv_2, hopf2d_mod_c_2 = continuation(hopf2d_modified,[-1 0],6,"beta",2:-0.05:-1,method = "pseudo_arclength",discretisation = "shooting",);

# Create traces
p1_u1 = scatter(x = algebraic_pv_1,y = algebraic_c_1[:,1],mode="lines",name="u1",showlegend=true, line=attr(color="blue"))
p1_u2 = scatter(x = algebraic_pv_1,y = algebraic_c_1[:,2],mode="lines",name="u2",showlegend=true, line=attr(color="red"))

p2_u1 = scatter(x = hopf2d_pv_1,y = hopf2d_c_1[:,1],mode="lines",name="u1",showlegend=false, line=attr(color="blue"))
p2_u2 = scatter(x = hopf2d_pv_1,y = hopf2d_c_1[:,2],mode="lines",name="u2",showlegend=false, line=attr(color="red"))

p3_u1 = scatter(x = hopf2d_mod_pv_1,y = hopf2d_mod_c_1[:,1],mode="lines",name="u1",showlegend=false, line=attr(color="blue"))
p3_u2 = scatter(x = hopf2d_mod_pv_1,y = hopf2d_mod_c_1[:,2],mode="lines",name="u2",showlegend=false, line=attr(color="red"))

p4_u1 = scatter(x = algebraic_pv_2,y = algebraic_c_2[:,1],mode="lines",name="u1",showlegend=false, line=attr(color="blue"))
p4_u2 = scatter(x = algebraic_pv_2,y = algebraic_c_2[:,2],mode="lines",name="u2",showlegend=false, line=attr(color="red"))

p5_u1 = scatter(x = hopf2d_pv_2,y = hopf2d_c_2[:,1],mode="lines",name="u1",showlegend=false, line=attr(color="blue"))
p5_u2 = scatter(x = hopf2d_pv_2,y = hopf2d_c_2[:,2],mode="lines",name="u2",showlegend=false, line=attr(color="red"))

p6_u1 = scatter(x = hopf2d_mod_pv_2,y = hopf2d_mod_c_2[:,1],mode="lines",name="u1",showlegend=false, line=attr(color="blue"))
p6_u2 = scatter(x = hopf2d_mod_pv_2,y = hopf2d_mod_c_2[:,2],mode="lines",name="u2",showlegend=false, line=attr(color="red"))


# Plot
p1 = plot([p1_u1, p1_u2], Layout(title="Algebraic NPC", xaxis_title = "c", yaxis_title = ""))
p2 = plot([p2_u1, p2_u2], Layout(title="Hopf NPC", xaxis_title = "beta", yaxis_title = ""))
p3 = plot([p3_u1, p3_u2], Layout(title="Modified Hopf NPC", xaxis_title = "beta", yaxis_title = ""))
p4 = plot([p4_u1, p4_u2], Layout(title="Algebraic PA", xaxis_title = "c", yaxis_title = ""))
p5 = plot([p5_u1, p5_u2], Layout(title="Hopf PA", xaxis_title = "beta", yaxis_title = ""))
p6 = plot([p6_u1, p6_u2], Layout(title="Modified Hopf PA", xaxis_title = "beta", yaxis_title = ""))


ncm_comp_plot = [p1 p2 p3; p4 p5 p6]
relayout!(ncm_comp_plot, height=500, width=800, title_text="Numerical Continuation Method Comparison")
ncm_comp_plot