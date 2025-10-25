function allJmax = JmaxTLSphasesT(A, phases, Edomain, tau_f, Nt, Nt_ts, Ncheb, tol)
    Nphases = length(phases);
    allJmax = zeros(1, Nphases);
    H0 = [0 0; 0 1];
    sigma_x = [0 1; 1 0];
    for iphases = 1:Nphases
        psi = SemiGlobalH(@(psi, t, v) H0*v - A*sin(t + phases(iphases))*sigma_x*v,...
            @(psi1, t1, psi2, t2) A*(sin(t2 + phases(iphases)) - sin(t1 + phases(iphases))).*(sigma_x*psi1), 1, [], Edomain, [1; 0], [0  tau_f], Nt, Nt_ts, Ncheb, tol);
        allJmax(iphases) = psi(2, end).*conj(psi(2, end));
    end
end