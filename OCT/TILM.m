function [field, psi, relE, conv, niter] = TILM(psi0, target, E, Epenal, T, dt, tol)
% OCT with a time independent Lagrange multiplier.
    Nt = T/dt;
%    lastfield = zeros(1, Nt);
    chi = target;
    psisol = ode45(@TILMpsi, [0 T], psi0, [], E, chi, Epenal);
    psi = deval(psisol, 0:dt:T);
%    lastfield = 1i/Epenal*(conj(chi(1))*psi(2, :) - conj(psi(1, :)).*chi(2));
    lastfield = 1/Epenal*imag(conj(chi(2))*psi(1, :) + conj(chi(1))*psi(2, :));
    relE = tol + 1;
    maxNiter = 1000;
    conv = zeros(1, maxNiter + 1);
    niter = 0;
    while relE>tol && niter<maxNiter
%         chi = (psi(:, Nt+1) + chi)/2;
%         chi = chi/norm(chi);
        overlap = target'*psi(:, Nt + 1);
        chi = overlap*target;
        niter = niter + 1;
        psisol = ode45(@TILMpsi, [0 T], psi0, [], E, chi, Epenal);
        psi = deval(psisol, 0:dt:T);
%        field = 1i/Epenal*(conj(chi(1))*psi(2, :) - conj(psi(1, :)).*chi(2));
        field = 1/Epenal*imag(conj(chi(2))*psi(1, :) + conj(chi(1))*psi(2, :));
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    niter
end