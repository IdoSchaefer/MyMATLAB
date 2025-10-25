function [Hoperations, fcouplingOp] = Hmats2Hops(H0, mu)
% Hop: A function handle of the form @(u, params, v, more arguments). Returns 
% H(u, params)v. params is a column vector which contains a list of
% time-dependent parameters evaluated at t.
% The "more arguments" will be written in the place of "varargin" in this program.
% Hdiff_op: A function handle of the form @(u1, params1, u2, params2, more arguments). Returns
% (H(u1, params1) - H(u2, params2))u1. u1 and params1 represent several time-points in 
% separate columns, while u2 and params2 represent a single time-point.
% Note: "more arguments" must be the same as for Hop.
   H0dag = H0';
   mu_dag = mu';
   Hoperations.psi = @(psi, field, v) H0*v - field*(mu*v);
   Hoperations.chi = @(chi, field, v) H0dag'*v - field*(mu_dag*v);
   Hoperations.diff_psi = @(psi1, field1, psi2, field2) (field2 - field1).*(mu*psi1);
   Hoperations.diff_chi = @(psi1, field1, psi2, field2) (field2 - field1).*(mu_dag*psi1);
   fcouplingOp = @(v) mu*v;
end