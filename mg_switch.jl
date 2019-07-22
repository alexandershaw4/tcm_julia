# voltage switch (magnesium block) for NMDA receptors
# Durstewitz, Seamans & Sejnowski 2000

function mg_switch(V)
    chv = 0.33*exp.(-0.06*V);
    s = 1.50265./(ones(1,8) + chv);
    return s
end
