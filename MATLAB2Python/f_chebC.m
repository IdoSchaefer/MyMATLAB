function Ccheb_f = f_chebC(t, Nz, Nt_ts, leftb, rightb, tol)
% The function computes the Chebyshev coefficients of \tilde{f}_{Nt_ts}(z, t),
% where z is the argument of the Chebyshev expansion, and t serves as a
% parameter.
% Input:
% t: Row vector of time-values
% Nz: The number of Chebyshev sampling points for the z argument
% Nt_ts: Defines the function as above
% leftb: The minimum boundary of the approximation domain
% rightb: The maximum boundary of the approximation domain
% tol: The tolerance parameter for the computation of \tilde{f}_{Nt_ts}(z, t)
% Output:
% Ccheb_f: A matrix containing the Chebyshev coefficients, where different
% t values are represented by separate columns.
    zsamp_cheb = cos(((1:Nz).'*2 - 1)*pi/(2*Nz));
    zsamp = 0.5*(zsamp_cheb*(rightb - leftb) + rightb + leftb);
    f_zt = f_fun(zsamp, t, Nt_ts, tol);
    Ccheb_f = chebcM(f_zt);
end
