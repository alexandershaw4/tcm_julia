# solve for a fixed point, such that fitted system responses (ERPs, SSRs) are
# expansions (deviations) about this point

function solvefixedpoint(P,M)

    M = copy(M);
    P = copy(P);

    ns = size(M['x'],1);
    np = size(M['x'],2);
    nk = size(M['x'],3);

    # void input
    M['u'] = [0];

    # initialise membrane potentials
    x = zeros(ns,np,nk);
    x[:,:,1] = repeat([-70],ns,np);
    x = float(copy(x));
    M['x'] = copy(x);

    a   = 2; # regularisation constant
    dnx = 0;

    for i = 1:128
        f,dfdx = tcm0(M['x'],M['u'],P,M);
        dfdx   = reshape(dfdx[:],ns*np*nk,ns*np*nk);
        dx     = -dfdx\f[:];

        ndx = norm(dx,Inf);
        if ndx < dnx
            a = a/2;
        end
        dnx = copy(ndx);

        # check convergence
        rdx = M['x'][:] + exp(-a)*dx[:];
        M['x'] = reshape(copy(rdx),ns,np,nk);

        if dnx < 1e-12;
            break
        end
    end

    x = copy(M['x']);
    return x

end
