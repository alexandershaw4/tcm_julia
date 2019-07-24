# julia implementation of thalamo-cortical model using Morris-Lecar-like
# conductance equations
#
# AS
import DataStructures
include("mg_switch.jl")   # the magnesium (voltage) switch
include("Ncdf.jl")        # j.d. williams CumulDistFun
include("dic_to_vect.jl") # [Un]Vectorising functions
include("jaco.jl")        # numerical differentiation routine

function tcm0(x,u,P,M)
    # a wrapper on the tcm code that does 1st order differentiation for
    # computation of the delay operator

    # state dimensions
    ns = size(M['x'],1);
    np = size(M['x'],2);
    nk = size(M['x'],3);

    # append Jacobian (J'*J)
    J,f = jaco(tcm,x,u,P,M,1);

    # append delay operator[s]
    D  = [.6 16];
    d  = -D.*exp.(P['D'])/1000;
    Sp = kron(ones(nk,nk),kron(Diagonal(ones(np)),Diagonal(ones(ns))));
    Ss = kron(ones(nk,nk),kron(ones(np,np),Diagonal(ones(ns))));

    # binary connectivity
    A = [1  0  1  0  0  1  0  1
         1  1  1  0  0  0  0  0
         1  1  1  0  0  0  0  0
         0  1  0  1  1  0  0  0
         0  0  0  1  1  0  0  0
         0  0  0  1  1  1  0  1
         0  0  0  0  0  0  1  1
         0  0  0  0  0  1  1  1];

    d0 = exp.(P['t']);
    Tc = zeros(np,np);
    Tc[7:8,1:6] = repeat([60*d0[1]],2,6);
    Tc[1:6,7:8] = repeat([20*d0[2]],6,2);
    Tc          = -Tc/1000;
    Tc          = Tc .* A;
    Tc = kron(ones(nk,nk),kron(Tc,Diagonal(ones(ns))));

    # concatenated inverted delay matrix
    Dp = (Ss.==0).*Ss;
    Ds = ((Sp.==0).*Sp) + (Ss.>0);
    D  = d[2]*Dp + d[1]*Ds + Tc;

    # Delay operator [Q]:
    #                 dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
    #                          = Q*f = Q*J*x(t)
    Q = inv( Diagonal(ones(56)) - D.*J );

    # Return order - these go for integration
    return f, J, Q
end

function tcm(x, u, P, M)
    # the main state equations for the conductance model

    # state space dimensions
    ns = size(M['x'],1);
    np = size(M['x'],2);
    nk = size(M['x'],3);
    x  = reshape(x,size(M['x']));

    # extrinsics, modulations & inputs
    A  = exp.(P['A']);
    AN = exp.(P['N'])
    B  = exp.(P['B']);
    BN = exp.(P['n'])
    C  = exp.(P['C']);

    # intrinsic connectivity
    G = exp.(P['H']);

    # connectivity switches
    #==========================================================================
    # 1 - excitatory spiny stellate cells (granular input cells)
    # 2 - superficial pyramidal cells     (forward  output cells)
    # 3 - inhibitory interneurons         (intrisic interneuons)
    # 4 - deep pyramidal cells            (backward output cells)
    # 5 - deep interneurons
    # 6 - thalamic projection pyramidal cells (with m- and h- currents)
    # 7 - thalamic reticular cells
    # 8 - thalamic relay cells (with m- and h- currents)
    =#

    # Extrinsic switches
    SA = [1 0 0 0 0
          0 1 0 0 0
          0 1 0 0 0
          1 0 0 0 0
          0 0 0 0 0
          0 0 0 0 0
          0 0 0 0 0
          0 0 0 0 0];
    SNMDA = SA;

    # Intrinsic: excitatory
    GEa = [0 0 0 0 0 2 0 2   # ss
           4 0 0 0 0 0 0 0   # sp
           4 4 0 0 0 0 0 0   # si
           0 4 0 0 0 0 0 0   # dp
           0 0 0 4 0 0 0 0   # di
           0 0 0 2 0 0 0 1/4 # tp
           0 0 0 0 0 0 0 2   # rt
           0 0 0 0 0 2 0 0]; # rl

    # Intrinsic: inhibitory
    GIa = [8  0  2  0  0  0  0  0  # ss
           0  16 16 0  0  0  0  0  # sp
           0  0  32 0  0  0  0  0  # si
           0  0  0  8  8  0  0  0  # dp
           0  0  0  0  16 0  0  0  # di
           0  0  0  0  8  8  0  0  # tp
           0  0  0  0  0  0  32 0  # rt
           0  0  0  0  0  0  8  32];# rl

    GEn = GEa; # NMDA   = same as AMPA
    GIb = GIa; # GABA-B = same as GABA-A

    # Receptor time-constants
    T  = P['T'];
    KE = exp(T[1])*1000/4;
    KI = exp(T[2])*1000/16;
    KN = exp(T[3])*1000/100;
    KB = exp(T[4])*1000/200;

    # Reversal potentials [voltages, mV]
    VL  = -70;
    VE  =  60;
    VI  = -90;
    VR  = -40;
    VN  =  10;
    VB  = -100;

    # Membrane capactiance(s)
    CV = P['V'];
    CV = diag(exp(Diagonal(CV)));
    CV = CV .* [128 32 32 128 64 128 256 64*8]';

    # K leak conductance
    GL = 1;

    # Firing approximation
    Vx = exp(P['S'][]) * 32;
    m  = Ncdf(x[:,:,1],VR,Vx);

    # Local extrinsic effects
    a  = zeros(ns,5);
    an = zeros(ns,5);

    a[:,1] = A*m[2];
    a[:,2] = A*m[4];
    a[:,3] = A*m[6];
    a[:,4] = A*m[7];
    a[:,5] = A*m[8];

    an[:,1] = AN*m[2];
    an[:,2] = AN*m[4];
    an[:,3] = AN*m[6];
    an[:,4] = AN*m[7];
    an[:,5] = AN*m[8];

    # Average background
    BE = exp.(P['B']) * 0.8

    # Inputs
    U = C.*u

    # new state vector
    f = float(copy(x));

    for i = 1:ns

        # multi-inputs
        dU = u[1]*( C[i]*[1, 1/64, 1/128, 1/128] );

        # dCurrents * firing
        E     = G[:,:,i]*GEa*m[i,:];
        ENMDA = G[:,:,i]*GEn*m[i,:];
        I     = G[:,:,i]*GIa*m[i,:];
        IB    = G[:,:,i]*GIb*m[i,:];

        # extrinsic coupling: excitatory
        E     = ( E     + repeat(BE,8,1) + SA*a[i,:] ) * 2;
        ENMDA = ( ENMDA + repeat(BE,8,1) + SNMDA*an[i,:] ) * 2;

        # exogenous inputs
        ic        = [8 1 2 4];
        E[ic]     = E[ic] + dU'
        ENMDA[ic] = ENMDA[ic] + dU'

        # Voltage equation
        f[i,:,1] = (GL.* repeat([float(VL)],8,1)-x[i,:,1]+
                x[i,:,2].*repeat([float(VE)],8,1)-x[i,:,1]+
                x[i,:,3].*repeat([float(VI)],8,1)-x[i,:,1]+
                x[i,:,5].*repeat([float(VB)],8,1)-x[i,:,1]+
                x[i,:,4].*repeat([float(VN)],8,1)-x[i,:,1].*
                        (mg_switch(x[i,:,1]'))')./CV;

        # Conductance equations
        f[i,:,2] = (E - x[i,:,2])*KE[i];
        f[i,:,3] = (E - x[i,:,3])*KI[i];
        f[i,:,4] = (E - x[i,:,4])*KN[i];
        f[i,:,5] = (E - x[i,:,5])*KB[i];

        DV = (1/[2,1,1,2.2,1,2,1,2])';
        DV = DV.*exp.(P['d']);

        f[i,:,2] = f[i,:,2].*DV;
        f[i,:,3] = f[i,:,3].*DV;
        f[i,:,4] = f[i,:,4].*DV;
        f[i,:,5] = f[i,:,5].*DV;

    end # end regions loop

    return f
end
