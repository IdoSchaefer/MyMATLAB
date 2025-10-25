function [field, psi, relE, conv, niter] = TLSTDtar(psi0, target, E, fguess, Epenal, tweight,T, dt, tol)
% fguess is a function handle of the form: @(t).
% tweight is a function handle of the form: @(t).
    Nt = T/dt;
%    field = zeros(1, Nt+1);
    lastfield = zeros(1, Nt+1);
%     chi = zeros(2, Nt+1);
%     psi = zeros(2, Nt+1);
    psisol = ode45(@TLSRK, [0 T], psi0, [], E, fguess);
%    psiT = deval(psisol, T);
    for ti = 1:(Nt+1)
        lastfield(ti) = fguess((ti-1)*dt);
    end
    if nargin>3
        field = lastfield;
        for ti = 1:(Nt + 1)
            wval = tweight((ti-1)*dt);
        end
        wval = wval/T;
    end
    relE = tol + 1;
%    conv = zeros(2, 1001);
    conv = zeros(1, 1001);
    maxNiter = 1000;
    niter = 0;
    while relE>tol && niter<maxNiter
%         overlap = target'*psiT;
%         chiT = overlap*target;
        if nargin>3
            psi = deval(psisol, 0:dt:T);
            overlap = sum(conj(target).*psi);
            conv(:, niter+1) = sum(wval.*overlap.*conj(overlap) - Epenal*field.*conj(field))*dt;
        end
        chisol = ode45(@TLSTDchi, [T 0], [0;0], [], E, psisol, Epenal, twheight);
        chi = deval(chisol, 0:dt:T);
%         for ti = 1:(Nt+1)
%             chi(:, ti) = RK4interp(chisol.y(:, end:-1:1), chisol.x(end:-1:1), (ti - 1)*dt);
%         end
        psisol = ode45(@TLSRKpsi, [0 T], psi0, [], E, chisol, Epenal);
        psi = deval(psisol, 0:dt:T);
%         for ti = 1:(Nt+1)
%             psi(:, ti) = RK4interp(psisol.y, psisol.x, (ti - 1)*dt);
%         end
%        psiT = psi(:, Nt + 1);
        field = 1i/Epenal*(conj(chi(1, :)).*psi(2, :) - conj(psi(1, :)).*chi(2, :));
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    if nargin>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = sum(conj(target).*psi);
        conv(:, niter+1) = sum(wval.*overlap.*conj(overlap) - Epenal*field.*conj(field))*dt;
        conv = conv(:, 1:(niter + 1));
    end
    niter
end