function [field, psi, relE, conv, niter] = TLSTDtarf(psi0, target, E, fguess, Epenal, T, dt, tol, tweight)
% target is a function handle of the form: @(t) [f(t);g(t)].
% fguess is a function handle of the form: @(t).
% tweight is a function handle of the form: @(t). The normalization over
% time to 1 is performed in the program.
    Nt = T/dt;
%    field = zeros(1, Nt+1);
    lastfield = zeros(1, Nt+1);
%     chi = zeros(2, Nt+1);
%     psi = zeros(2, Nt+1);
%    chiT = zeros(2, 1);
    tarT = target(T);
    psisol = ode45(@TLSRK, [0 T], psi0, [], E, fguess);
    psiT = deval(psisol, T);
    for ti = 1:(Nt+1)
        lastfield(ti) = fguess((ti-1)*dt);
    end
    if nargin > 8        
        % normalization to 1:
        integw = quad(tweight, 0, T);
        tweight = @(t) tweight(t)/integw;
        wval = zeros(1, Nt + 1);
        for ti = 1:(Nt + 1)
            wval(ti) = tweight((ti-1)*dt);
        end
    else
        tweight = @(t) 1/T;
        wval = ones(1, Nt + 1)/T;
    end
%        wval = wval/T;
    maxNiter = 1000;
    if nargin>3
        field = lastfield;
        conv = zeros(1, maxNiter + 1);
        tarval = zeros(2, Nt + 1);
        for ti = 1:(Nt + 1)
            tarval(:, ti) = target((ti-1)*dt);
        end        
    end
    relE = tol + 1;
%    conv = zeros(2, 1001);
    niter = 0;
    while relE>tol && niter<maxNiter
         overlap = tarT'*psiT;
         chiT = overlap*tarT;
%        chiT = [0; 1i*lastfield(Nt + 1)/(conj(psiT(1))*Epenal)];
%         chiT(1) = -1i*conj(lastfield(Nt + 1))/Epenal*(tarT(1) - psiT(1))/(conj(tarT(2))*psiT(1) - tarT(1)*conj(psiT(2)));
%         chiT(2) = (conj(chiT(1))*psiT(2) + 1i*lastfield(Nt + 1)/Epenal)/conj(psiT(2));
        if nargin>3
            psi = deval(psisol, 0:dt:T);
%            overlap = sum(conj(target).*psi);
%             overlap = sum(conj(tarval).*psi);
%             conv(:, niter+1) = sum(wval.*overlap.*conj(overlap) - Epenal*field.*conj(field))*dt;
            overlap_allt = sum(conj(tarval).*psi);
            conv(:, niter+1) = overlap*conj(overlap) + sum(wval.*overlap_allt.*conj(overlap_allt) - Epenal*field.*conj(field))*dt;
        end
        chisol = ode45(@TLSTDchif, [T 0], chiT, [], E, psisol, Epenal, target, tweight);
        chi = deval(chisol, 0:dt:T);
%         for ti = 1:(Nt+1)
%             chi(:, ti) = RK4interp(chisol.y(:, end:-1:1), chisol.x(end:-1:1), (ti - 1)*dt);
%         end
        psisol = ode45(@TLSRKpsi, [0 T], psi0, [], E, chisol, Epenal);
        psi = deval(psisol, 0:dt:T);
%         for ti = 1:(Nt+1)
%             psi(:, ti) = RK4interp(psisol.y, psisol.x, (ti - 1)*dt);
%         end
        psiT = psi(:, Nt + 1);
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
%        overlap = sum(conj(target).*psi);
%         overlap = sum(conj(tarval).*psi);
%         conv(:, niter+1) = sum(wval.*overlap.*conj(overlap) - Epenal*field.*conj(field))*dt;
        overlap = tarT'*psiT;
        overlap_allt = sum(conj(tarval).*psi);
        conv(:, niter+1) = overlap*conj(overlap) + sum(wval.*overlap_allt.*conj(overlap_allt) - Epenal*field.*conj(field))*dt;
        conv = conv(:, 1:(niter + 1));
    end
    niter
end