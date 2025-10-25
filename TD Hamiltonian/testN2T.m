function [result, relE] = testN2T(f, xdomain, N, testpoints)
% The program tests the procudure of conversion between Newton and Taylor
% polynomial coefficients, in order to test its numerical stability for
% high orders.
% Input:
% f: a function handle. represents the function to be interpolated.
% xdomain: The approximation domain
% N: The degree of the approximation
% testpoints: The points in which the result is computed
% Output:
% result: The resulting approximation in all testpoints
% relE: The relative error in all testpoints
    xcheb = -cos((0:N)*pi/N);
    xdlength = xdomain(2) - xdomain(1);
    x = (xcheb*xdlength + xdomain(1) + xdomain(2))/2;
    Cr2t = r2Taylor4(x, xdlength);
    Cnewton = devdif(x*4/xdlength, f(x));
    Ctaylor = zeros(1, N + 1);
    Ctaylor(1) = Cnewton(1);
    for Newtoni = 2:(N + 1)
        Ctaylor(1:Newtoni) = Ctaylor(1:Newtoni) + Cnewton(Newtoni)*Cr2t(Newtoni, 1:Newtoni);
    end
    result = zeros(size(testpoints));
    for polydeg = 0:N
        result = result + Ctaylor(polydeg + 1)*testpoints.^polydeg;
    end
    relE = (result - f(testpoints))./f(testpoints);
end