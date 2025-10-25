function [allNt, allmv, aller, allCPU] =  error_decayRK7(fderiv, t0tf, u0, uf_exact, minNt, Nsamp, varargin)
% The function computes the data of the error decay with reduction of the
% time-step for RK7.
% Input:
% fderiv: a function handle of the form @(t, u, varargin). The value of the function
% is the derivative of u with respect to t.
% t0tf: a vector that contains initial and final time, in the form: [t0 tf]
% u0: The initial values of the solution vector u.
% uf_exact: The exact solution at the final time (obtained analytically or
% numerically by a highly accurate propagation)
% Output:
% allNt: The number of equally spaced time-points in the propagation
% process for all sampling points in the error decay curve
% allmv: The number of matrix-vector multiplications (or Hamiltonian
% operations) for all sampling points in the error decay curve
% aller: The relative errors from the provided uf_exact for all sampling
% points in the error decay curve
    allNt = zeros(Nsamp, 1);
    allmv = zeros(Nsamp, 1);
    aller = zeros(Nsamp, 1);
    allCPU = zeros(Nsamp, 1);
    t_interval_length = (t0tf(2) - t0tf(1));
    for degi = 1:Nsamp
        deg = log10(minNt) + (degi-1)*0.1;
        Nt = round(10^deg);
        allNt(degi) = Nt;
        dt = t_interval_length/Nt;
        tic
        u = RK7uf(fderiv, t0tf, u0, dt, varargin{:});
        allCPU(degi) = toc;
        matvecs = Nt*9;
        allmv(degi) = matvecs;
        if ~isfinite(u)
            fprintf('\nError.\n')
        end
        error = norm(u - uf_exact)/norm(uf_exact);
        aller(degi) = error;
    end
    figure
    plot(log10(allmv), log10(aller), '-o')
    xlabel('log(matvecs)')
    ylabel('log(error)')
end