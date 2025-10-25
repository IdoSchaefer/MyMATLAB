function [Phi_out, all_norm_res, alldiagPhi] = smooth(updateM, L, omega, h, Phi_in, fv, Niter)
    all_norm_res = ones(1, Niter + 1);
    N = sqrt(length(Phi_in)) + 1;
    if nargout>2
        alldiagPhi = zeros(N - 1, Niter + 1);
        alldiagPhi(:, 1) = Phi_in(1:N:(N - 1)^2);
    end
    Phi_out = Phi_in;
    res = fv - L*Phi_out;
    all_norm_res(1) = sqrt(sum(res.^2))/(N - 1);
    for iteri = 1:Niter
        Phi_out = updateM*Phi_out - h^2*fv/(2*(2 + omega));
        res = fv - L*Phi_out;
        all_norm_res(iteri + 1) = sqrt(sum(res.^2))/(N - 1);
        if nargout>2
            alldiagPhi(:, iteri + 1) = Phi_out(1:N:(N - 1)^2);
        end
    end
end