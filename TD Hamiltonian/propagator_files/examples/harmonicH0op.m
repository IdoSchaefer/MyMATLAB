function result = harmonicH0op(v, x, K)
% The function returns the harmonic oscillator H0 operation on a vector v (m=1).
% x: The x grid.
% K: The kinetic energy vector in the p domain. Has to be computed in the
% following way:
%     p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
%     p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
%     K = p.^2/2;
    result = ifft(K.*fft(v)) + x.^2/2.*v;
end