function dctM = dctIM(M)
% The function performes a discrete cosine transform of the first kind.
% Suitable if it is desired to include the boundaries of the domain in
% the transform.
% The convention used for the factor of the sums, makes the transform its
% own inverse.
    dim = size(M);
    N = dim(1);
    dctM = 1/sqrt(2*(N - 1))*fft([M; M((N - 1):-1:2, :)]);
    dctM = dctM(1:N, :);
    if isreal(M)
        dctM = real(dctM);
    elseif isreal(1i*M)
        dctM = 1i*imag(dctM);
    end    
end