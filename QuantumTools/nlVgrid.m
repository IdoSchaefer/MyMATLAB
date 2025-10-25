function [xdomain, dx, x] = nlVgrid(Vf, N, m)
% The function computes the optimal Fourier grid of N points, for a
% non-linear potential. The grid satisfies the following demands:
% 1. The potential value is equal at the boundaries of the grid (roughly,
%   because the grid doesn't include the right boundary).
% 2. The maximal potential energy is equal to the maximal kinetic energy,
%   P^2/(2*m), where Pmax = pi/dx.
% The 2 demands determine the potential uniquely, for a given number of
% grid points, N.
% input:
% Vf: a function handle of x, represents the potential function. Vf(0) = 0 
% m: the mass. The default is: m = 1.
% output:
% xdomain: a 2 element vector, [xmin xmax].
% x: the grid
    if nargin < 3
        m = 1;
    end
    % The above mentioned demands impose the following demand on xmax:
    xmaxfun =  @(x) Vf(x) - Vf(x - N*pi/sqrt(2*m*Vf(x)));
    xmax = fzero(xmaxfun, 10);
    % and the following demand on dx:
    dx = pi/sqrt(2*m*Vf(xmax));
    xmin = xmax - N*dx;
    xdomain = [xmin, xmax];
    x = xmin:dx:(xmax - dx);
end