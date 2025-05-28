function [rho,u] = rhoNu(f,ksi)
    % This function is self explainitory

    rho = sum(f, 1);

    u = (pagemtimes(ksi,f))./rho;
end

