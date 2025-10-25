function Dv = DcosinesMap(v, r)
% The function returns the derivative of a function sampled at the 
% points r, using mapped cosine differentiation.
% v: A vector that contains the function values. It is assumed to satisfy the
% cosine series boundary conditions (zero derivative).
% r: A vector that contains the sampling points. The mapping function r(x) should satisfy the
% cosine series boundary conditions (zero derivative).
    dim = size(v);
    if dim(1) == 1
        N = dim(2) - 1;
        v = v.';
        r = r.';
    else
        N = dim(1) - 1;
    end
    Dv = zeros(N + 1, 1);
    k = (0:N).';
    dctv = dctI(v);
    dctr = dctI(r);
    Dv(2:N) = dstI(-k(2:N).*dctv(2:N))./dstI(-k(2:N).*dctr(2:N));
    if dim(1) == 1
        Dv = Dv.';
    end
end
