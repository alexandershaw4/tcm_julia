# Compute the 1st or 2nd order partial (numerical) derivates of a function
#
# usage: [j,ip] = jaco(fun,x0,V,verbose,order)
#
# (order 1:) Compute the 1st order partial derivatives (gradient)
# of a function using:
#   j(ip,:) = ( f(x(ip)+h)  - f(x(ip)-h) )  / (2 * h)
#
# (order 2:) Compute the 2nd order derivatives (curvature):
#
# j(ip,:) = [ (f0 - f1) / 2 / d ] ./ [ (f0 - 2 * fx + f1) / d ^ 2 ]
#
#
# if order==1, when j is square, it is the Jacobian
# if order==2, when j is square, it is the Hessian
#

import DataStructures

function(fun,x,u,P,M,order)
    f0 = eval(fun(x,u,P,M));
    fx = f0;
    j  = zeros( length(x), length(x) );
    V  = f0;

    for i = 1:length(x)

        x0 = vec(x);
        x1 = vec(x);
        d  = x[i] + V[i]*exp(-8);

        if d == 0
            d = 0.01;
        end

        x0[i] = x0[i] + d;
        x1[i] = x1[i] - d;

        f0 = eval(fun,reshape(x0,size(x)),u,P,M);
        f1 = eval(fun,reshape(x1,size(x)),u,P,M);
        j[i,:] = (f0 - f1) / (2d);

        if order == 2
            deriv1 = (f0 - f1) / 2 / d;
            deriv2 = (f0 - 2 * fx + f1) / d^2;
            j[i,:] = feriv1 ./ deriv2;
        end
    end

    j = j';

    return j
end
