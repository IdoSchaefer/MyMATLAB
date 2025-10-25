function [Gop, Gdiff_op] = Hmats2Gop_dipole(H0, mu, f_field)
% The function constructs the input function handles for the Hamiltonian
% operations for the program SemiGlobal.m, for a Hamiltonian with a time-dependent
% dipole interaction term of the following form:
% H(u, t) = H0 - epsilon(u, t)mu
% The program applies to the case that H0 and mu are given as 
% time-independent matrices, and epsilon(u, t) is given as the functional 
% form of the field.
% Input:
% H0: The drift Hamiltonian matrix
% mu: The dipole operator matrix
% f_field: A function handle of the form: @(u, t); returns the field at
% several time points.
%     f_field input: u: The time-dependent state at the time points specified by t;
%                       separate time point are represented by separate columns.
%                    t: A row vector of the time points for evaluation of
%                    the field
%     f_field output: A row vector of the field evaluated at the time-points
%                     specified by t.
% Output:
% Gop: A function handle of the form @(u, t, v). Returns  G(u, t)v,
% where G(u, t) = -(1i/hbar)H(u, t). In the current program we
% use the units: hbar = 1.
% Gdiff_op: A function handle of the form @(u1, t1, u2, t2). Returns
% (G(u1, t1) - G(u2, t2))u1. u1 and t1 represent several time-points in 
% separated columns, while u2 and t2 represent a single time-point.
   Gop = @(u, t, v) -1i*(H0*v - f_field(u, t)*(mu*v));
   Gdiff_op = @(u1, t1, u2, t2) 1i*(f_field(u1, t1) - f_field(u2, t2)).*(mu*u1);
end