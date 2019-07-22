
function dict_to_vec(d)
    v = float;
    for k in keys(d)
        v = vcat(v, vec(float(d[k])))
    end
    return v
end

function vec_to_dict(Px,d0)
    NP = copy(P)

    if Px[1] == float
        # remove useless first entry
        Px = Px[2:end];
    end

    for k in keys(d0)

        # get size of this key/field/entry
        S0 = size(P[k]);
        Sn = prod(S0);

        # copy it into place
        NP[k] = reshape( Px[1:Sn], S0 );

        # remove these values from beginning of vector
        Px = Px[Sn+1:end];
    end
    return NP
end
