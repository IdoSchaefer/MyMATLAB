function dctv = dctI0b(v)
% The function computes the dctI of v, which contains the dctI coefficients
% for 0 boundaries. The two first coefficients are not included in v, and
% are determined by the boundary constraints.
    dim = size(v);
    if dim(1) == 1
        N = dim(2) + 2;
        v = v.';
    else
        N = dim(1) + 2;
    end
    % v_div_d is v./d, where d represents the vector of the d's in a dctI sum,
    % with d_0 = d_(N-1) = 2, and all the rest are 1.
    v_div_d = v;
    v_div_d(N - 2) = 0.5*v_div_d(N - 2);
    v1 = -sum(v_div_d(2:2:(N - 2)));
    dctv = 1/sqrt(2*(N - 1))*fft([-2*sum(v_div_d(1:2:(N - 2))); v1; v; v((N - 3):-1:1); v1]);
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