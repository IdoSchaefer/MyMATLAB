function [allNt, allmv, aller, all_est_ers] = error_forced_harmonic1(Nt_ts, Nfm, minNt, Nsamp, Niter, Arnoldi)
    if nargin<5
        Niter = 1;
    end
    if nargin<6
        Arnoldi = false;
    end
    % Constructing the grid:
    L = 16*sqrt(pi);
    Nx = 128;
    dx = L/Nx;
    x = (-L/2:dx:(L/2 - dx)).';
    % Constructing the kinetic energy matrix diagonal in the p domain:
    p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    K = p.^2/2;
    % The potential energy matrix diagonal in the x domain:
    V = x.^2/2;
    % The harmonic oscillator ground state:
    fi0 = pi^(-1/4)*exp(-x.^2/2)*sqrt(dx);
    T = 10;
    % Computing the expectation value of x in all the time points:
    mx_ex = (-0.5*sin(T)*T);
    % Computing the expectation value of p in all the time points:
    mp_ex = -0.5*(sin(T) + T*cos(T));
    angle_analytic = (T/2 - (sin(2*T)/4 - T*cos(2*T)/2)/8);
    Uex_phase = pi^(-1/4)*exp(-1i*angle_analytic)*exp(1i*(mp_ex.*(x - mx_ex/2)) - (x - mx_ex).^2/2)*sqrt(dx);
    if Arnoldi
        [allNt, allmv, aller, all_est_ers] = error_decaySG2(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v),...
            @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [], fi0, Uex_phase, [0 T], Nt_ts, Nfm, minNt, Nsamp, Niter);
    else
        [allNt, allmv, aller, all_est_ers] = error_decaySG2(@(u, t, v) -1i*Hpsi(K, V + x*cos(t), v),...
            @(u1, t1, u2, t2) -1i*x*(cos(t1) - cos(t2)).*u1, 0, [], [-195*1i, 0], fi0, Uex_phase, [0 T], Nt_ts, Nfm, minNt, Nsamp, Niter);
    end
end
    
