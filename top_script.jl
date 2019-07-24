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

# copy in some parameter values
P['C'] = [-0.1374 -0.0046 -0.0019 -0.0013];
P['V'] = [-0.0383 0.0022 0.0201 0.0128 0.0186 -0.0370 7.8399e-04 9.1260e-04]';
P['T'] = [0.3111 -0.0585 0.6164 3.3609];
P['L'] = [0.2571];
P['d'] = [-0.0042 -0.6699 0.4010 0.0066 -0.0067 0.1795 -0.1221 -0.0767]';
P['t'] = [-0.0043 0.0012]';
P['H'][:,:,1] =
[0.3525         0    0.0100         0         0         0         0    0.2173
-0.0368   -0.0060    0.1421         0         0         0         0         0
     0   -9.4951   -0.0581         0         0         0         0         0
     0         0         0   -8.7795         0         0         0         0
     0         0         0    0.1224   -0.1636         0         0         0
     0         0         0   -0.0298         0    0.0753         0    0.0006
     0         0         0         0         0         0    0.0153         0
     0         0         0         0         0    0.0412         0   -0.1671];
P['A'] = [0];
P['N'] = [0];
P['B'] = [0];
P['n'] = [0];

# solve for fixed point [solvefixedpoint.jl]
x = solvefixedpoint(P,M);

M['x'] = x;

# check the integration & transfer functions [integrate.jl]
k, w, y_out, y, pst = integrate(P,M);

# fit the model using AO curvature optimisation
