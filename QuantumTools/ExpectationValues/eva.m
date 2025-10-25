function ma = eva(psi, a0, field)
% The function computes the expectation value of the acceleration of psi in
% all time steps. The system is under the influence of a time dependent
% field. The dipole approximation in the x gauge is used.
% psi: The state in all time steps.
% a0: A vector which represents the acceleration form of the time independent Hamiltoniam. Assumed to
% be a function of x.
% field: The external field in all time steps. The default is zero field.
    Nt = size(psi, 2);
    if nargin<3
        field = zeros(1, Nt);
    end
    ma = zeros(1, Nt);
    for ti = 1:Nt
        ma(ti) = psi(:, ti)'*((a0 + field(ti)).*psi(:, ti));
    end
end