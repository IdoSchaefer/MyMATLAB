function [field, psi, relE, conv, niter, mniterc] = OCchebMLS(psi0, target, E0, Edomain, miu, fguess, Epenal, T, dt, Nt_ts, Ncheb, tol)
% fguess is a function handle of the form: @(t).
% E0 is a vector containing the eigenenergies of the time independent
% Hamiltonian matrix.
    Nt = T/dt;
    dim = length(psi0);
    H0 = diag(E0);
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    t_ts = 0.5*(tcheb+1)*dt;
%    field = zeros(1, Nt+1);
%    lastfield = zeros(1, Nt+1);
%    chi = zeros(dim, Nt+1);
    psi = zeros(dim, Nt+1);
    allmniter = 0;
tic
%    [allpsi, lastfield, mniter] = solveOCMLS(@VtguessMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt);
    [allpsi, lastfield, mniter] = solveOCMLSih(@ihguessMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, fguess, t_ts, dt);
    allmniter = allmniter + mniter;
toc
    psiT = allpsi(:, Nt_ts, Nt);
%     for ti = 1:(Nt+1)
%         lastfield(ti) = fguess((ti-1)*dt);
%     end
    if nargin>3
        field = lastfield;
    end
%    lastfield = field;
    relE = tol + 1;
    maxNiter = 1000;
    conv = zeros(1, maxNiter + 1);
    niter = 0;
    while relE>tol && niter<maxNiter
        tic
        overlap = target'*psiT;
        chiT = overlap*target;
        if nargin>3
            conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        end
%        [allchi, field1, mniter] = solveOCMLS(@VtchiMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-2, allpsi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allchi, field1, mniter] = solveOCMLSih(@ihchiMLS, H0, Edomain, miu, chiT, [T 0], Nt, Nt_ts, Ncheb, tol*1e-2, allpsi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        allmniter = allmniter + mniter;
%        [allpsi, field, mniter] = solveOCMLS(@VtpsiMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        [allpsi, field, mniter] = solveOCMLSih(@ihpsiMLS, H0, Edomain, miu, psi0, [0 T], Nt, Nt_ts, Ncheb, tol*1e-2, allchi(:, Nt_ts:-1:1, Nt:-1:1), Epenal);
        allmniter = allmniter + mniter;
        psiT = allpsi(:, Nt_ts, Nt);
        relE = norm(field - lastfield)/norm(field);
        lastfield = field;
        niter = niter + 1;
        toc
    end
    if niter==maxNiter
        fprintf('The program has failed to achieve the desired tolerance.\n')
    end
    psi(:, 1:Nt) = allpsi(:, 1, :);
    psi(:, Nt + 1) = allpsi(:, Nt_ts, Nt);
    if nargin>3        
%       conv(:, niter+1) = psi(:, Nt + 1).*conj(psi(:, Nt + 1));
        overlap = target'*psiT;
        conv(:, niter+1) = overlap*conj(overlap) - Epenal*sum(field.^2)*dt;
        conv = conv(:, 1:(niter + 1));
    end
    niter
    mniterc = allmniter/(2*niter + 1)
end