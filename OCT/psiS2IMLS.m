function psiI = psiS2IMLS(psiS, E0, T)
% The function converts the wave functions in all time points, from the
% Schrodinger picture (S) to the interaction picture (I). Intended to the
% energy representation.
% psiS: the wave function in the Schrodinger picture in all time points.
% E0: a vector containing the eigenenergies of H0
% T: the final time
    Npsi = size(psiS, 1);
    Nt = size(psiS, 2) - 1;
    dt = T/Nt;
    % if E0 is a row vector, then transpose:
    if size(E0, 2) ~= 1
        E0 = E0.';
    end
    psiI = zeros(Npsi, Nt + 1);
    for ti = 1:(Nt + 1)
        psiI(:, ti) = exp(1i*E0*(ti - 1)*dt).*psiS(:, ti);
    end
end