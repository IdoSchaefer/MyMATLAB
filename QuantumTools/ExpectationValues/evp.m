function exval_p = evp(psix, p)
% The function computes the time-dependent expectation value of the
% momentum exval_mu from the time-dependent state vector in the x domain 
% psix and the momentum vector in the p domain, p.
% psi: a matrix; the state vector in the x domain at all time points, where
% different time-points are represented by separate columns.
% p: a column vector; the momentum vector, represented in the p domain.
% in the x domain
    exval_p = sum(conj(psix).*ifft(p.*fft(psix)));
%     Nt = size(psi, 2);
%     mp = zeros(1, Nt);
%     for ti = 1:Nt
%         mp(ti) = psi(:, ti)'*ifft(p.*fft(psi(:, ti)));
%     end
end