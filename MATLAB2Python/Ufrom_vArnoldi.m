function [U, f_error, fUerror] = Ufrom_vArnoldi(v_vecs, timeM, Upsilon, RvKr, samplingp, capacity, Nt_ts, Nfm, tol)
    % The function computes the solution from the v vectors.
    % Input:
    % timeM: The matrix of the t powers for the required time-points
    % Output:
    % U: The solution in the required time points
    % f_error: The estimated relative error of the computation of f_fun(G,t(Nt_ts))v_vecs(:, Nkr + 1)
    % fUerror: The estimated relative error of U resulting from the
    % computational error of f(\tilde{G},t(Nt_ts))v_vecs(:, Nkr + 1)
        
    % The required time-points are given by their first power:
    tp = timeM(2, :);
    % Computing the divided differences for the Newton expansion of f_fun, for
    % all time-points:
    fdvd_ts = divdif(samplingp.'/capacity, f_fun(samplingp, tp, Nt_ts, tol).').';
    % The computation of the Newton expansion of f_fun(G, tp)v_vecs(:, Nkr + 1)
    % in the Krylov space, for all tp:
    fGv_kr = RvKr(1:Nfm, :)*fdvd_ts;
    % Computing the solution in all tp: 
    U = v_vecs(:, 1:(end-1))*timeM + Upsilon(:, 1:Nfm)*fGv_kr;        
    if nargout>1
        % The absolute error:
        f_error_abs = abs(fdvd_ts(Nfm + 1, Nt_ts - 1))*norm(RvKr(:, Nfm + 1));
        % The relative error:
        f_error = f_error_abs/norm(fGv_kr(:, Nt_ts - 1));
        % The relative error of U, resulting from this computation:
        fUerror = f_error_abs/norm(U(:, Nt_ts - 1));
    end
end
