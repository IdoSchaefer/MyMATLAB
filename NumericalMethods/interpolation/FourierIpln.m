function vnew = FourierIpln(v, xdomain, xnew)
% The function computes an interpolation of the equally sampled signal v
% in the points xnew. The signal is interpolated using the Fourier basis
% functions. xnew needn't be equally spaced. v should satisfy the Fourier
% periodic boundary conditions.
% xdomain: The boundaries of the Fourier domain of v.
    N = length(v);
    L = xdomain(2) - xdomain(1);
    k = 0:2*pi/L:2*pi*(N - 1)/L;
    k(ceil((N+1)/2):N) = k(ceil((N+1)/2):N) - 2*pi*N/L;
    rows = size(v, 1);
    if rows == 1
        v = v.';
    end
    if size(xnew, 1) == 1
        xnew = xnew.';
    end
    vk = fft(v);
    vnew = (exp(1i*(xnew - xdomain(1))*k)*vk)/N;
    if rows == 1
        vnew = vnew.';
    end
end