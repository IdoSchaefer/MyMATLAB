function [leftb, rightb, ev] = Re_evdomain(M)
% The function finds the real eigenvalue domain [leftb rightb] of a matrix
% M.
% ev: The eigen
    [P, D] = eig(M);
    ev = diag(D);
    Re_ev = real(ev);
    [Re_ev, ev_indices] = sort(Re_ev);
    leftb = min(Re_ev);
    rightb = max(Re_ev);
    ev = ev(ev_indices);
end