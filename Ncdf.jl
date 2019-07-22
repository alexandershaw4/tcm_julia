# JD Williams approximation

function Ncdf(x,u,v)
    x = (x .- u)/sqrt(abs(v));
    F = sqrt.(1 .- exp.(-(2/pi)*x.^2))/2;
    i = x .< 0;
    F[i] = -F[i];
    F = F .+ 1/2;
    return F
end
