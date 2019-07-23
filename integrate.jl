# Numerical integration and spectral response from a
# conductance based model [transfer function]

using FFTW;
using DSP;
using Dierckx;

include("tcm.jl");
include("jaco.jl");
include("solvefixedpoint.jl");
include("Setup_tcm.jl");
include("dic_to_vect.jl");
include("integrate.jl");
include("Ncdf.jl");
include("mg_switch.jl");

function integrate(P,M)

    # settings
    w   = 4:80;
    x   = M['x'];
    dt  = 1/1200;
    Fs  = 1/dt;
    tn  = 2;
    pst = (0:dt:tn-dt)*1000;

    # input (over time)
    mu = exp.(P['R']);
    drive = repeat(mu,1,length(pst));

    # solve for a fixed point
    x = solvefixedpoint(P,M);
    M['x'] = copy(x);

    # vectorise x
    x = x[:];

    # initialise full states time series
    y = zeros(ns*np*nk,length(pst));

    # precompute delay operators
    fx, dfdx, D = tcm0(M['x'],0,P,M);

    #
    n       = ns*np*nk;
    N       = 2;
    Q       = (expm1.(dt*D*dfdx/N) - Diagonal(ones(n))).*inv(dfdx);

    # start point
    v = copy(fx[:]);

    # Runge-Kutta (1st order) numerical integration w/ Delays
    for i = 1:length(pst)
        for j = 1:2
            global v += real(Q) * tcm(v,drive[i],P,M)[:] ;
        end
        # record as expansion about fp
        y[:,i] = v - x;
    end

    # Compute forward projection: J'*y
    LFP = exp.(P['J'][:])'*y;

    # Use FFTW FFT
    y0 = rfft(y') |> fftshift;
    fq = rfftfreq(length(pst), 1200) |> fftshift

    # 1d spline to freqs of interest
    Contributing = exp.(P['J'][:]) .!= 0;
    Contributing = findall(Contributing);
    y_out        = zeros(length(Contributing),length(w));

    for i = 1:length(Contributing)
        y_out[i,:] = Spline1D(fq,real(y0[:,Contributing[i]]),k=1)(w);
    end

    k = sum(y_out,dims=1);

    return k w y_out y pst
end
