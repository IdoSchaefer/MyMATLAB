function [fieldw, psi, evmiut, evmiuw, relE, conv, niter, mniterc] = OCfMLSK(psi0, E0, Edomain, miu, fguess, filterE, filtermiu, T, dt, Nt_ts, Ncheb, tol)
% fguess: a function handle. A function of the time.
% E0 is a vector containing the eigenenergies of the time independent
% Hamiltonian matrix.
    Nt = T/dt;
    dw = pi/T;
    Npsi = length(psi0);
    H0 = diag(E0);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
    fieldt = zeros(1, Nt+1);
    maxNiter = 100;
    conv = zeros(1, maxNiter + 1);
%    lastfield = zeros(1, Nt+1);
%    chi = zeros(Npsi, Nt+1);
    psi = zeros(Npsi, Nt + 1);
    chiT = zeros(Npsi, 1);
    vfilterE = zeros(1, Nt+1);
    vfiltermiu = zeros(1, Nt);
    for wi = 1:Nt
        vfilterE(wi) = filterE((wi-1)*dw);
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
%        fieldw(wi) = fguess((wi-1)*dw);
    end
    vfilterE(Nt + 1) = filterE(Nt*dw);
%    fieldt = dct(fieldw);
    allmniter = 0;
tic
    [allpsi, fieldt(1:Nt), mniter] = solveOCMLSih(@ihguessMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt);
toc
    allmniter = allmniter + mniter;
    fieldt(Nt + 1) = fguess(T);
    lastfieldw = dctI(fieldt);
    [chiihterm, evmiuw] = get_chiihterm(allpsi, miu, vfiltermiu);
%    [evmiuw, evmiufil] = get_evmiuf(allpsi, miu, vfiltermiu);
    if nargin>5
        conv(1) = Jval(evmiuw(1:Nt), vfiltermiu(1:Nt), lastfieldw(1:Nt), vfilterE(1:Nt), dw);
    end
tic
    [allchi, fieldt1, mniter] = solveOCMLSih(@ihchigMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt, chiihterm);
toc
    allmniter = allmniter + mniter;
%    chimiupsi = get_cmp(allchi(:, Nt_ts:-1:1, Nt:-1:1), allpsi, miu);
    chimiupsi = get_cmp(allchi(:, Nt_ts:-1:1, Nt:-1:1), allpsi, miu);
tic
    [allpsi, chimiupsi, mniter] = solveOCfMLS(@ihpsifMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, ...
            chimiupsi, allchi(:, Nt_ts:-1:1, Nt:-1:1), vfilterE);
toc
    allmniter = allmniter + mniter;
    fieldw = dctI(chimiupsi).*vfilterE;
    relE = norm(fieldw - lastfieldw)/norm(fieldw);
    lastfieldw = fieldw;
    [chiihterm, evmiuw] = get_chiihterm(allpsi, miu, vfiltermiu);
%    [evmiuw, evmiufil] = get_evmiuf(allpsi, miu, vfiltermiu);
    if nargin>5
        conv(2) = Jval(evmiuw(1:Nt), vfiltermiu(1:Nt), fieldw(1:Nt), vfilterE(1:Nt), dw);
    end
%    relE = tol + 1;
    niter = 1;
    while relE>tol && niter<maxNiter
        tic
        [allchi, chimiupsi, mniter] = solveOCfMLS(@ihchifMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-2, ...
            chimiupsi, allpsi(:, Nt_ts:-1:1, Nt:-1:1), vfilterE, chiihterm);
        allmniter = allmniter + mniter;
        [allpsi, chimiupsi, mniter] = solveOCfMLS(@ihpsifMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, ...
            chimiupsi, allchi(:, Nt_ts:-1:1, Nt:-1:1), vfilterE);
        allmniter = allmniter + mniter;
        fieldw = dctI(chimiupsi).*vfilterE;
        relE = norm(fieldw - lastfieldw)/norm(fieldw)
        lastfieldw = fieldw;
        [chiihterm, evmiuw] = get_chiihterm(allpsi, miu, vfiltermiu);
%        [evmiuw, evmiufil] = get_evmiuf(allpsi, miu, vfiltermiu);
        niter = niter + 1;
        if nargin>5
            conv(niter + 1) = Jval(evmiuw(1:Nt), vfiltermiu(1:Nt), fieldw(1:Nt), vfilterE(1:Nt), dw);
        end
        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    evmiut = evmiuMLS(psi, miu);
    niter
    conv = conv(1:niter+1);
    mniterc = allmniter/(2*niter + 1)
end

function [chiihterm, evmiuw] = get_chiihterm(allpsi, miu, vfiltermiu)
    [Npsi, Nt_ts, Nt] = size(allpsi);
    chiihterm = zeros(Npsi, Nt_ts, Nt);
    [evmiuw, evmiufil] = get_evmiuf(allpsi, miu, vfiltermiu);
    for tsi = 1:Nt
        for ti = 1:Nt_ts
            chiihterm(:, ti, tsi) = -evmiufil(Nt_ts - ti + 1, Nt - tsi + 1)*miu*allpsi(:, Nt_ts - ti + 1, Nt - tsi + 1);
        end
    end
end

function [evmiuw, evmiufil] = get_evmiuf(allpsi, miu, vfiltermiu)
    [Npsi, Nt_ts, Nt] = size(allpsi);
    evmiufil = zeros(Nt_ts, Nt);
    allevmiut = zeros(Nt_ts-1, Nt + 1);
    for tsi = 1:Nt
        for ti = 1:Nt_ts-1
            allevmiut(ti, tsi) = allpsi(:, ti, tsi)'*miu*allpsi(:, ti, tsi);
        end
    end
    allevmiut(1, Nt + 1) = allpsi(:, Nt_ts, Nt)'*miu*allpsi(:, Nt_ts, Nt);
    for ti = 1:Nt_ts-1
% Note that we use here a different dct from the one used for the field.
% The reason is, that the number of points we use here is Nt, so the
% boundary isn't included.
        evmiuw = dct(allevmiut(ti, 1:Nt));
        evmiufil(ti, :) = idct(evmiuw.*vfiltermiu);
    end
    evmiuw = dct(allevmiut(1, 2:(Nt+1)));
    evmiufil(Nt_ts, :) =  idct(evmiuw.*vfiltermiu);
% For analysis purposes, it may be important to use a dct which represents
% correctly the time grid, and self consistent with the dct used for the
% field:
    evmiuw = dctI(allevmiut(1, :));
end

function J = Jval(evmiuw, vfiltermiu, fieldw, vfilterE, dw)
% The function computes the value of the functional J.
    notsmall = (abs(fieldw)>1e5*eps) + (abs(vfilterE)>1e5*eps);
    J2fun = fieldw(notsmall == 2).^2./vfilterE(notsmall == 2);
%    J2fun = J2fun(isfinite(J2fun));
%    J2fun = J2fun(abs(J2fun)<100);
%    conv(:, niter) = (sum(0.5*evmiuw.^2.*vfiltermiu) - sum(fieldw(vfilterE~=0).^2./vfilterE(vfilterE~=0)))*dw;
    J = (sum(0.5*evmiuw.^2.*vfiltermiu) - sum(J2fun))*dw;
end

function chimiupsi = get_cmp(chi, psi, miu)
    [Npsi, Nt_ts, Nt] = size(psi);
    chimiupsi = zeros(1, Nt + 1);
%     for tsi = 1:Nt
%         for ti = 1:Nt_ts
%             chimiupsi(ti, tsi) = -imag(allchi(:, ti, tsi)'*miu*allpsi(:, ti, tsi));
%         end
%     end
     for tsi = 1:Nt
         chimiupsi(tsi) = -imag(chi(:, 1, tsi)'*miu*psi(:, 1, tsi));
     end
     chimiupsi(Nt+1) = -imag(chi(:, Nt_ts, Nt)'*miu*psi(:, Nt_ts, Nt));
end