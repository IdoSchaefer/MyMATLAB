function dctv = dctI0lb(v)
% The function computes the dctI of v, which contains the dctI coefficients
% for a 0 left boundary. The last coefficient is not included in v, and is
% determined by the boundary constraint.
    dim = size(v);
    if dim(1) == 1
        N = dim(2) + 1;
        v = v.';
    else
        N = dim(1) + 1;
    end
    dctv = 1/sqrt(2*(N - 1))*fft([v; -2*(v(1)/2 + sum(v(2:(N - 1)))); v((N - 1):-1:2)]);
    dctv = dctv(1:N);
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end    
end