# Set up model parameters and data structures
using DataStructures
using LinearAlgebra

function Init_Params()

    ns = 1; # num regions
    np = 8; # num cells per region
    nk = 7; # num states per cell

    M = OrderedDict{Char,Array{Float64}}()
    P = OrderedDict{Char,Array{Float64}}()

    # Fill out M struct
    M['x'] = zeros(1,8,7);

    # Fill out Parameter struct
    P['A'] = ones(ns,ns); # Extrns connectivity [ampa]
    P['N'] = ones(ns,ns); # Extrns connectivity [nmda]
    P['B'] = ones(ns,ns); # trial-spec extrns connectivity [ampa]
    P['n'] = ones(ns,ns); # trial-spec extrns connectivity [nmda]
    P['C'] = ones(ns,1);  # inputs
    P['D'] = [0 0];       # delays [intr], [extr]
    P['E'] = [0];         # background activity
    P['G'] = [0];         # trial spec intrinsics (off)
    P['H'] = zeros(np,np,ns);

    # State contributions to output
    P['J'] = ones(ns,np,nk)*-1000;
    P['J'][1,[1 2 4 6],1] = (log(Diagonal([.2,.8, .2, .2]))).diag';
    P['L'] = [4];

    # Other parameters
    P['S']  = [0];              # firing
    P['T']  = zeros(ns,4);      # decay rates: AMPA, GABAA, NMDA, GABAB
    P['V'] = [0,0,0,0,0,0,0,0]; # membrane capacitance
    P['d'] = [0,0,0,0,0,0,0,0]; # conductance delays
    P['t'] = [0,0];             # T->C & C->T delays
    P['R'] = [0];               # d.c. input

    return P, M
end
