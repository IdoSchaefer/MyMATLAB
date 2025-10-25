function [field, psi, relE, conv, niter] = OCTLSK(psi0, target, E, fguess, Epenal, T, dt, tol)
% fguess is a function handle of the form: @(t).
    Nt = T/dt;
%    field = zeros(1, Nt+1);
%    chi = zeros(2, Nt+1);
    fguess = @(u, t) fguess(t);
    [psi, lastfield] = TLSexpH(E, fguess, psi0, T, Nt);
    if nargin>3
        field = lastfield;
    end
%     fchi = @(chit, t) 1i/Epenal*(conj(chit(1))*psi(2, Nt - t/(-dt) + 1) - conj(psi(1, Nt - t/(-dt) + 1))*chit(2));
%     fpsi = @(psit, t) 1i/Epenal*(conj(chi(1, Nt - t/dt + 1))*psit(2) - conj(psit(1))*chi(2, Nt - t/dt + 1));
%    lastfield = field;
    relE = tol + 1;
%    conv = zeros(2, 1001);
    conv = zeros(1, 1001);
    maxNiter = 1000;
    niter = 0;
    while relE>tol && niter<maxNiter
%        chiT = (target'*psi(:, Nt + 1))*target;
        overlap = target'*psi(:, Nt + 1);
        chiT = overlap*target;
        if nargin>3
%            conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));            
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.*conj(field))*dt;
        end
%        chi = TLSexpH(E, fchi, chiT, -T, Nt);
        chi = TLSKchi(E, chiT, psi, Epenal, T, Nt);
%        [psi, field] = TLSexpH(E, fpsi, psi0, T, Nt);
        [psi, field] = TLSKpsi(E, psi0, chi, Epenal, T, Nt);
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
%   conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
    overlap = target'*psi(:, Nt + 1);
    conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.*conj(field))*dt;
    conv = conv(:, 1:(niter + 1));
    niter
end