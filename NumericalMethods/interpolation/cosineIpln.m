function vnew = cosineIpln(v, xdomain, xnew)
% The function computes an interpolation of the equally sampled signal v
% in the points xnew. The signal is interpolated using cosine basis
% functions. xnew needn't be equally spaced. v should satisfy the cosine
% series boundary conditions of zero derivative.
% xdomain: The boundaries of the domain of v.
    N = length(v) - 1;
    L = xdomain(2) - xdomain(1);
    k = 0:pi/L:N*pi/L;
    rows = size(v, 1);
    if rows == 1
        v = v.';
    end
    if size(xnew, 1) == 1
        xnew = xnew.';
    end
    vk = dctI(v);
    vnew = sqrt(2/N)*(cos((xnew - xdomain(1))*k)*[0.5*vk(1); vk(2:N); 0.5*vk(N + 1)]);
    if rows == 1
        vnew = vnew.';
    end
end