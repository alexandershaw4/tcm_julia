# top-level script for tcm.
#
# actions:

# requirements
include("tcm.jl");
include("jaco.jl");
include("solvefixedpoint.jl");
include("Setup_tcm.jl");
include("dic_to_vect.jl");
include("integrate.jl");
include("Ncdf.jl");
include("mg_switch.jl");

# collect some input data

# set up parameters [Setup_tcm.jl]

# solve for fixed point [solvefixedpoint.jl]

# check the integration & transfer functions [integrate.jl]

# fit the model using AO curvature optimisation
