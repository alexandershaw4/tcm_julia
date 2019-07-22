# julia implementation of thalamo-cortical model using Morris-Lecar-like
# conductance equations
#
# AS
import DataStructures

function tcm(x::float64, u::float64, P, M)

    # state space dimensions
    ns = size(M['x'],1)
    np = size(M['x'],2)
    nk = size(M['x'],3)
    x  = reshape(x,size(M['x']))

    # extrinsics, modulations & inputs
    A = exp(P['A'])
    B = exp(P['B'])
    C = exp(P['C'])

    # intrinsic connectivity
    G = exp(P['H'])

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
          0 0 0 0 0]
    SNMDA = SA

    # Intrinsic: excitatory
    GEa = [0 0 0 0 0 2 0 2   # ss
           4 0 0 0 0 0 0 0   # sp
           4 4 0 0 0 0 0 0   # si
           0 4 0 0 0 0 0 0   # dp
           0 0 0 4 0 0 0 0   # di
           0 0 0 2 0 0 0 1/4 # tp
           0 0 0 0 0 0 0 2   # rt
           0 0 0 0 0 2 0 0]  # rl

    # Intrinsic: inhibitory
    GIa = [8  0  2  0  0  0  0  0  # ss
           0  16 16 0  0  0  0  0  # sp
           0  0  32 0  0  0  0  0  # si
           0  0  0  8  8  0  0  0  # dp
           0  0  0  0  16 0  0  0  # di
           0  0  0  0  8  8  0  0  # tp
           0  0  0  0  0  0  32 0  # rt
           0  0  0  0  0  0  8  32]# rl

    GEn = GEa # NMDA   = same as AMPA
    GIb = GIa # GABA-B = same as GABA-A

    # Receptor time-constants
    T  = P['T']
    KE = exp(T[1])*1000/4
    KI = exp(T[2])*1000/16
    KN = exp(T[3])*1000/100
    KB = exp(T[4])*1000/200

    # Reversal potentials [voltages, mV]
    VL  = -70
    VE  =  60
    VI  = -90
    VR  = -40
    VN  =  10
    VB  = -100

    # Membrane capactiance(s)
    CV = P['V']
    CV = diag(exp(Diagonal(CV)))
    CV = CV .* [128 32 32 128 64 128 256 64*8]'

    # K leak conductance
    GL = 1

    # Firing approximation
    Vx = exp(P['S'][]) * 32
    m  = Ncdf(x[:,:,1],VR,Vx)

    # Local extrinsic effects
    a  = zeros(ns,5)
    an = zeros(ns,5)







end
