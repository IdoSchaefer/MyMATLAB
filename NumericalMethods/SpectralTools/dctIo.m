function dctv = dctIo(v)
% The function performes a discrete cosine transform of the first kind.
% Suitable if it is desired to include the boundaries of the domain in
% the transform.
% The convention used for the factor of the sums, makes the transform its
% own inverse.
% The convention for the factors of the 1'st and last terms, makes the
% tranformation orthogonal.
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
    else
        N = dim(1);
    end
    dctv = 1/sqrt(2*(N - 1))*fft([v(1)*sqrt(2); v(2:N-1); v(N)*sqrt(2) ; v((N - 1):-1:2)]);
    dctv = dctv(1:N);
    dctv([1, N]) = dctv([1, N])/sqrt(2);
    if dim(1) == 1
        dctv = dctv.';
    end
end