function [U, fUerror] = Ufrom_vCheb(v_tpols, timeM, Vcheb, Ccheb_f, f_error)
% The function computes the solution for the Chebyshev algorithm at all
% time points specified by the transpose Vandermonde matrix timeM.
% Input:
% v_tpols: The v_j vectors excluding the last one, j=0,1,...,Nt_ts-1, in
% seperate columns; the vector coefficients of the Taylor time-polynomials
% in the solution equation
% timeM: The matrix of the t powers for the required time-points
% Vcheb: The T_k(\tilde{G})v_{Nt_ts} vectors, k=0,1,...,Nt_ts-1, in
% seperate columns
% Ccheb_f: The Chebyshev coefficients of \tilde{f}_{Nt_ts}(z, t) in the
% required time-points specified by timeM, as computed by the function f_chebC
% Output:
% U: The computed solution at the required time-point
% fUerror: Estimation of the error resulting from the Chebyshev approximation
% for the computation of the function of matrix
    fGv = Vcheb*Ccheb_f;
    U = v_tpols*timeM + fGv;
    if nargout>1
        fUerror = f_error*norm(fGv(:, end - 1))/norm(U(:, end - 1));
    end
end
