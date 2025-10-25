function [T, R, total_norm, normL_m, normR_p] = TRcoef(E, Vf, xdomain, Nx)
% The function computes the transmission and reflection coefficients for a
% potential in a given spatial domain.
% Input:
% E: The energy of the system
% Vf: A function handle; the potential function. It is possible to insert
% the potential value vector itself instead.
% xdomain: The domain of the x grid
% Nx: The number of eqidistant points in the x grid. Needed to define the
% accuracy of the computation.
% Output:
% T: The transmission coefficient
% R: The reflection coefficient
% total_norm: The total norm of the reflected and transmitted waves. Should
% be 1 for a real potential.
    dx = (xdomain(2) - xdomain(1))/(Nx - 1);
    x = xdomain(1):dx:xdomain(2);
    if length(Vf) == 1
        % If Vf is a function handle:
        V = Vf(x);
    else
        % If Vf is a vector:
        V = Vf;
    end
    if E == V(1)
        E = V(1) + 10*eps;
    end
    k_i = sqrt(2*(E - V(1)));
    k_L = k_i;
    % Computation of the transfer matrix M:
    M = eye(2);
    for xi = 2:Nx
        if E~=V(xi)
            k_R = sqrt(2*(E - V(xi)));
        else
            k_R = sqrt(2*(10*eps));
            %k_R = 1e5*eps;
        end
        M = ([(k_L + k_R)*exp(1i*(k_L - k_R)*x(xi))   (k_R - k_L)*exp(-1i*(k_L + k_R)*x(xi));
              (k_R - k_L)*exp(1i*(k_L + k_R)*x(xi))   (k_R + k_L)*exp(-1i*(k_L - k_R)*x(xi))]*M/(2*k_R));
        k_L = k_R;        
    end
    L_minus = -M(2, 1)/M(2, 2);
    R_plus = det(M)/M(2, 2);
    % det(M) is a difference. When using complex potential, the M values
    % may be very large, and it will result in a numerical instability.
    % Maybe TRcoef1 performes better.
    normL_m = L_minus*conj(L_minus);
    normR_p = R_plus*conj(R_plus);
    T = R_plus*conj(R_plus)*real(k_R)/real(k_i);
    R = L_minus*conj(L_minus);
    total_norm = T + R;
end