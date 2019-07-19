# julia implementation of thalamo-cortical model using
# Morris-Lecar-like conductance equations
# AS
import DataStructures

function tcm(x::float64, u::float64, P, M)

    ns = size(M['x'],1)
    np = size(M['x'],2)
    nk = size(M['x'],3)
