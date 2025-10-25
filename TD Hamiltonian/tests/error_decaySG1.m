function [allNt, allmv, aller, est_ers, allCPU] = error_decaySG1(Gop, Gdiff_op, Gdiff_matvecs, ihfun, ev_domain, ui, uf_exact,...
    t_interval, Nt_ts, Nfm, minNt, Nsamp, Niter, varargin)
% The function computes the data of the error decay with reduction of the
% time-step for the semi-global propagator. The procedure SemiGlobalDirty.m is
% tentatively used (this should be fixed after completion of a reliable error
% estimation code). The number of propagation time-points in the output data increases
% nearly linearly in a logaritmic scale.
% Input:
% Gop: A function handle of the form @(u, t, v, more arguments). Returns 
% G(u, t)v (in quantum mechanics, G(u, t) = -iH(u, t)/hbar). The "more arguments"
% will be written in the place of "varargin" in this program.
% Gdiff_op: A function handle of the form @(u1, t1, u2, t2, more arguments). Returns
% (G(u1, t1) - G(u2, t2))u1. u1 and t1 represent several time-points in 
% separated columns, while u2 and t2 represent a single time-point.
% Note: "more arguments" must be the same as for Gop.
% Gdiff_matvecs: The number of matrix-vector multiplications required for
% the computation of the operation of Gdiff_op for each time-point (usually less than 2).
% (The number of matrix-vector multiplications is counted as the number of
% large scale operations with the highest scaling with the dimension of the
% problem.)
% ihfun: A function handle of the form: @(t, more arguments). Returns the
% inhomogeneous term vector. If there is no inhomogeneous term, put []
% instead.
% Note: "more arguments" must be the same as for Gop.
% ev_domain: The boundaries of the eigenvalue domain of G, when a Chebyshev
% algorithm is used for the computation of the function of matrix. For an
% Arnoldi algorithm, put [].
% ui: the initial state vector.
% uf_exact: The exact solution at the final time (obtained analytically or
% numerically by a highly accurate propagation)
% t_interval: The time interval of the propagation problem; vector of the form:
% [ti, tf]
% Nt_ts: the number of interior Chebyshev time points in each time-step, used during
% the computational process (M in the article).
% Nfm: the number of expansion terms for the computation of the function of
% matrix (K in the article).
% minNt: The minimal number of time-points in the error decay curve
% Nsamp: The number of sampling points in the error decay curve
% Niter: The fixed number of iterations in each time-step, excluding the
% first, which is iterated to the machine precision
% Output:
% allNt: The number of equally spaced time-points in the propagation
% process for all sampling points in the error decay curve
% allmv: The number of matrix-vector multiplications (or Hamiltonian
% operations) for all sampling points in the error decay curve
% aller: The relative errors from the provided uf_exact for all sampling
% points in the error decay curve
% max_ers: The maximal estimated errors from the different error sources
% (the max_errors output structure of SemiGlobal.m) for all sampling points
% in the error decay curve
    if nargin<13
        Niter = 1;
    end
    allNt = zeros(Nsamp, 1);
    allmv = zeros(Nsamp, 1);
    aller = zeros(Nsamp, 1);
    allCPU = zeros(Nsamp, 1);
    est_ers.texp_cheap = zeros(Nsamp, 1);
    est_ers.texp_exact = zeros(Nsamp, 1);
    est_ers.fU = zeros(Nsamp, 1);
    est_ers.conv_cheap = zeros(Nsamp, 1);
    est_ers.conv_exact = zeros(Nsamp, 1);
    est_ers.total = zeros(Nsamp, 1);
    for degi = 1:Nsamp
        % Computing the number of time-points Nt for a particular sampling
        % point in the error decay curve; note that log10(Nt) increases
        % approximately linearly with degi.
        deg = log10(minNt) + (degi-1)*0.1;
        Nt = round(10^deg);
        allNt(degi) = Nt;
        tic
        [u, ~, matvecs, ~, history] = SemiGlobalDirty(Gop, Gdiff_op, Gdiff_matvecs,...
            ihfun, ev_domain, ui, t_interval, Nt, Nt_ts, Nfm, eps, Niter, 16, false, varargin{:});
        allCPU(degi) = toc;
        allmv(degi) = matvecs;
        est_ers.texp_cheap(degi) = history.total_errors.texp_cheap;
        est_ers.texp_exact(degi) = history.total_errors.texp_exact;
        est_ers.fU(degi) = history.total_errors.fU;
        est_ers.conv_cheap(degi) = history.total_errors.conv_cheap;
        est_ers.conv_exact(degi) = history.total_errors.conv_exact;
        est_ers.total(degi) = history.total_errors.total;
        if ~isfinite(u(:,end))
            fprintf('\nError.\n')
        end
        error = norm(u(:, end) - uf_exact)/norm(uf_exact);
        aller(degi) = error;
    end
    figure
    plot(log10(allmv), log10(aller), '-o')
    xlabel('log(matvecs)')
    ylabel('log(error)')
end