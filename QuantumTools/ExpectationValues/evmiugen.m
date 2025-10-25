function mmiu = evmiugen(psi, x, q, w)
% The function computes the x expectation value in all time points.
% w: The weight function for the grid.
    Nt = size(psi, 2);
    mmiu = zeros(1, Nt);
    for ti = 1:Nt
        mmiu(ti) = psi(:, ti)'*(w.*x.*psi(:, ti));
    end
    if nargin>2
        mmiu = mmiu*q;
    end
end