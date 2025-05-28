function [f_eq] = eqm_d2q9(rho, u, ksi, w)

    f_eq = rho .* (1 + pagemtimes(ksi',u)/(1/3) + (pagemtimes(ksi',u).^2)/(2/9) - sum(u.^2, 1)/(2/3)) .* w';
    
end

