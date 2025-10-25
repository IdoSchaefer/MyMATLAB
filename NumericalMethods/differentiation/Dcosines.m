function Dv = Dcosines(v, domainL)
% The function returns the derivative of a function that its derivative at
% the boundaries is 0, and can be represented by a cosine series.
% v: A vector that contains the function values.
% domainL: The length of the domain.
    dim = size(v);
    if dim(1) == 1
        N = dim(2) - 1;
        v = v.';
    else
        N = dim(1) - 1;
    end
    Dv = zeros(N + 1, 1);
    k = (0:N).';
    dctv = dctI(v);
    Dv(2:N) = dstI(-k(2:N).*dctv(2:N))*pi/domainL;
    if dim(1) == 1
        Dv = Dv.';
    end
end