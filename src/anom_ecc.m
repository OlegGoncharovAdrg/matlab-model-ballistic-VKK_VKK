function [E] = anom_ecc(M,e)
    format long
    E = M;
    k = 1;
    err = 1e-10;
    while (k>err)
        y = E-(E-e*sin(E)-M)/(1-e*cos(E));
        k = abs(abs(E)-abs(y));
        E = y;
    end
end