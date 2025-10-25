function dstv = dstI(v)
% The function performes a discrete sine transform of the first kind.
% The convention used for the factor of the sums, makes the transform its
% own inverse.
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
    else
        N = dim(1);
    end
    dstv = 1i/sqrt(2*(N + 1))*fft([0; v; 0; -v(N:-1:1)]);
    dstv = dstv(2:(N + 1));
    if dim(1) == 1
        dstv = dstv.';
    end
    if isreal(v)
        dstv = real(dstv);
    elseif isreal(1i*v)
        dstv = 1i*imag(dstv);
    end
end