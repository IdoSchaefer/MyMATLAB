function [fieldw, psi, evmiut, evmiuw, relE, conv, niter] = OCfMLS(psi0, H0, miumat, fguess, filterE, filtermiu, T, dt, tol)
% fguess: a function handle. A function of the frequency.
    Nt = T/dt;
    dw = pi/T;
    Npsi = length(psi0);
%    W = pi/dt;
    fieldw = zeros(1, Nt + 1);
%    fieldt = zeros(1, Nt + 1);
    chimiupsi = zeros(1, Nt + 1);
%    lastfieldw = zeros(1, Nt + 1);
    chi = zeros(Npsi, Nt + 1);
%    chi0 = zeros(Npsi, 1);
    chiT = zeros(Npsi, 1);
%    chi0 = psi0;
    psi = zeros(Npsi, Nt + 1);
    vfilterE = zeros(1, Nt + 1);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
        fieldw(wi) = fguess((wi-1)*dw);
    end
    fieldt = dctI(fieldw);
    lastfieldw = fieldw;        
% tic
% %     psisol = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
% %    [t, psi] = ode45(@OCguess, 0:dt:T, psi0, [], K, Vf, x, fguess);
%     [tpsi, psisol] = ode45(@OCguess, [0 T], psi0, [], K, Vf, x, fguess);
%     psisol = psisol.';
% toc
%     for ti = 1:(Nt+1)
%         lastfield(ti) = fguess((ti-1)*dt);
%     end
%     if nargin>3
%         field = lastfield;
%     end
    relE = tol + 1;
    conv = zeros(1, 1001);
    maxNiter = 20;
    niter = 0;
    while relE>tol && niter<maxNiter
%        tic
        [tpsi, psisol] = ode45(@OCfMLSpsi, [0 T], psi0, [], H0, miumat, fieldt, Nt, dt);
        psisol = psisol.';
        for ti = 1:(Nt+1)
            psi(:, ti) = RK4interp(psisol, tpsi, (ti - 1)*dt);
        end        
        evmiut = evmiuMLS(psi, miumat);
        evmiuw = dctI(evmiut);
        evmiufil = dctI(evmiuw.*vfiltermiu);
%        chiT = -evmiut(Nt + 1)*miumat*psi(:, Nt + 1);
%        chiT = -H0*psi(:, Nt + 1);
%        psi = deval(psisol, 0:dt:T);        
%        chisol = ode45(@OCchi, [T 0], chiT, [], K, Vf, x, psisol, Epenal);
%        [t, chi] = ode45(@OCchi, T:-dt:0, chiT, [], K, Vf, x, psi, Epenal, T, dt);
%        [tchi, chisol] = ode45(@OCfMLSchi, [0 T], chi0, [], H0, miumat, fieldt, psi, evmiufil, Npsi, Nt, dt);
        [tchi, chisol] = ode45(@OCfMLSchi, [T 0], chiT, [], H0, miumat, fieldt, psi, evmiufil, Npsi, Nt, dt);
        chisol = chisol.';
        chisol = chisol(:, end:-1:1);
        tchi = tchi(end:-1:1);
        for ti = 1:(Nt+1)
            chi(:, ti) = RK4interp(chisol, tchi, (ti - 1)*dt);
%             chi(:, ti) = RK4interp(chisol.y(:, end:-1:1), chisol.x(end:-1:1), (ti - 1)*dt);
        end
%        chi = deval(chisol, 0:dt:T);
        for ti = 1:Nt
%             chit = RK4interp(chisol, tchi, (ti - 1)*dt);
%             psit = RK4interp(psisol, tpsi, (ti - 1)*dt);
            chimiupsi(ti) = -imag(chi(:, ti)'*miumat*psi(:, ti));
        end
%        fieldw = dct(chimiupsi).*vfilterE;
        fieldw = 0.95*fieldw + 0.05*dctI(chimiupsi).*vfilterE;
        fieldt = dctI(fieldw);
%         for ti = 1:(Nt+1)
%             field(ti) = -imag(chi(:, ti)'*(x.*psi(:, ti)))/Epenal;
%         end
        relE = norm(fieldw - lastfieldw)/norm(fieldw);
        lastfieldw = fieldw;
        niter = niter + 1;
        if nargin>5
            evmiut = evmiuMLS(psi, miumat);
            evmiuw = dctI(evmiut);
%            J2fun = fieldw.^2./vfilterE;
            notsmall = (abs(fieldw)>1e5*eps) + (abs(vfilterE)>1e5*eps);
            J2fun = fieldw(notsmall == 2).^2./vfilterE(notsmall == 2);
%            J2fun = J2fun(isfinite(J2fun));
            J2fun = J2fun(abs(J2fun)<100);
%            conv(:, niter) = (sum(0.5*evmiuw.^2.*vfiltermiu) - sum(fieldw(vfilterE~=0).^2./vfilterE(vfilterE~=0)))*dw;
            conv(:, niter) = (sum(0.5*evmiuw.^2.*vfiltermiu) - sum(J2fun))*dw;
        end
%        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
%     [tpsi, psisol] = ode45(@OCfpsi, [0 T], psi0, [], K, Vf, x, fieldt);
%     psisol = psisol.';
%     for ti = 1:(Nt+1)
%         psi(:, ti) = RK4interp(psisol, tpsi, (ti - 1)*dt);
%     end
%     evmiut = evmiu(psi, x);
%     evmiuw = dct(evmiut);
%    psi = deval(psisol, 0:dt:T);
    if nargin>5       
% %       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
%         overlap = target'*psiT;
%         conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        conv = conv(:, 1:niter);
    end
    niter
end