# top-level script for tcm.
#
# actions:

# requirements
using FFTW;
using DSP;
using Dierckx;
using DataStructures;
using LinearAlgebra;

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
P,M = Init_Params();

# solve for fixed point [solvefixedpoint.jl]
x = solvefixedpoint(P,M);

# check the integration & transfer functions [integrate.jl]
k, w, y_out, y, pst = integrate(P,M);

# fit the model using AO curvature optimisation
