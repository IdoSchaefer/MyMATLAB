function [allNt, allmv, aller, all_est_ers, allCPU] = error_decaySG2(Gop, Gdiff_op, Gdiff_matvecs, ihfun, ev_domain, ui, uf_exact,...
    t_interval, Nt_ts, Nfm, minNt, Nsamp, Niter, options, varargin)
% The function computes the data of the error decay with reduction of the
% time-step for the semi-global propagator. The procedure SemiGlobal1.m is
% used. The number of propagation time-points in the output data increases
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
% first, which is iterated Niter+3 times.
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
    if nargin<14        
        options = SGdefault_op;
        options.Niter = Niter;
        options.tol1st = eps;
        options.save_memory = true;
    end
    data = SGdata(options);
    allNt = zeros(Nsamp, 1);
    allmv = zeros(Nsamp, 1);
    aller = zeros(Nsamp, 1);
    allCPU = zeros(Nsamp, 1);
    all_est_ers.texp_cheap = zeros(Nsamp, 1);
    Nt_ts_is_even = mod(Nt_ts, 2) == 0;
    if ~Nt_ts_is_even
        all_est_ers.texp_cheap_odd = zeros(Nsamp, 1);
    end
    if Nt_ts_is_even || options.texp_er_odd
        all_est_ers.texp_exact = zeros(Nsamp, 1);
    end
    all_est_ers.fm = zeros(Nsamp, 1);
    all_est_ers.conv = zeros(Nsamp, 1);
    all_est_ers.conv_texp = zeros(Nsamp, 1);
    all_est_ers.conv_fm = zeros(Nsamp, 1);
    if options.conv_er_cheb
        all_est_ers.conv_cheb = zeros(Nsamp, 1);
    end
    all_est_ers.total = zeros(Nsamp, 1);
    for degi = 1:Nsamp
        % Computing the number of time-points Nt for a particular sampling
        % point in the error decay curve; note that log10(Nt) increases
        % approximately linearly with degi.
        deg = log10(minNt) + (degi-1)*0.1;
        Nt = round(10^deg);
        allNt(degi) = Nt;
        tic
        [u, ~, matvecs, est_errors] = SemiGlobal1(Gop, Gdiff_op, Gdiff_matvecs,...
            ihfun, ev_domain, ui, t_interval, Nt, Nt_ts, Nfm, [], options, data, varargin{:});
        allCPU(degi) = toc;
        allmv(degi) = matvecs.propagation;
        all_est_ers.texp_cheap(degi) = est_errors.texp_cheap;
        if ~Nt_ts_is_even
            all_est_ers.texp_cheap_odd(degi) = est_errors.texp_cheap_odd;
        end
        if Nt_ts_is_even || options.texp_er_odd
            all_est_ers.texp_exact(degi) = est_errors.texp_exact;
        end
        all_est_ers.fm(degi) = est_errors.fm;
        all_est_ers.conv(degi) = est_errors.conv;
        all_est_ers.conv_texp(degi) = est_errors.conv_texp;
        all_est_ers.conv_fm(degi) = est_errors.conv_fm;
        if options.conv_er_cheb
            all_est_ers.conv_cheb(degi) = est_errors.conv_cheb;
        end
        all_est_ers.total(degi) = est_errors.total;
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