function [T, check] = trans_coef(E, DV, Da)
    Njunctions = length(DV);
    k_i = sqrt(2*E);
    k_L = sqrt(2*E);
    % Computation of the transfer matrix M:
    V_R = 0;
    a = 0;
    M = eye(2);
    for junctioni = 1:Njunctions
        V_R = V_R + DV(junctioni);
        a = a + Da(junctioni);
        if E~=V_R
            k_R = sqrt(2*(E - V_R));
        else
            k_R = 10^5*eps;
        end
        M = ([(k_L + k_R)*exp(1i*(k_L - k_R)*a)   (k_R - k_L)*exp(-1i*(k_L + k_R)*a);
              (k_R - k_L)*exp(1i*(k_L + k_R)*a)   (k_R + k_L)*exp(-1i*(k_L - k_R)*a)]*M/(2*k_R));
        k_L = k_R;        
    end
    L_minus = -M(2, 1)/M(2, 2);
    R_plus = det(M)/M(2, 2);
    % Check:
    check = L_minus*conj(L_minus) + R_plus*conj(R_plus)*k_R/k_i - 1;
    T = R_plus*conj(R_plus)*k_R/k_i;
end