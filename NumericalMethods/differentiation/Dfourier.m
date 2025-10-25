function Dv = Dfourier(v, domainL)
% The function returns the derivative of a periodic function, using Fourier
% differentiation.
% v: A vector that contains the function values.
% domainL: The length of the domain.
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
    else
        N = dim(1);
    end
    dx = domainL/N;
    k = (0:(2*pi/domainL):(2*pi*(1/dx - 1/domainL))).';
    k((N/2 + 1):N) = k((N/2 + 1):N) - 2*pi/dx;
    Dv = ifft(1i*k.*fft(v));
    if dim(1) == 1
        Dv = Dv.';
    end
end