theta = 0:pi/16:pi;
% Sampling Chebyshev points in the domain [0, 2*pi]:
x = (-cos(theta) + 1)*pi;
% A vector function valued at the sampling points:
fx = cos((0.01:0.01:10).'*x);
% For the divdif_morag code:
xt = x.'; fxt = fx.';
% Regular iterative code:
tic, for j= 1:1000, a1 = divdif(x,fx); end, toc
% Morag's code, based on the solution of a lower triangular matrix:
tic, for j= 1:1000, a2 = divdif_morag(xt,fxt); end, toc
% A code based on a linear transformation of the function values into the
% coefficients:
tic, for j= 1:1000, a3 = divdif_mat(x,fx); end, toc
% Run this code several times, to see the typical times for each algorithm.


