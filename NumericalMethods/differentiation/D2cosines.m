function D2v = D2cosines(v, domainL)
% The function returns the 2'nd derivative of a function that its derivative at
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
    k = (0:N).';
    D2v = dctI(-k.^2.*dctI(v))*pi/domainL;
    if dim(1) == 1
        D2v = D2v.';
    end
end