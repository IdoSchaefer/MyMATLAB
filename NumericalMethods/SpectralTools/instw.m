function w = instw(signal, dx)
% The function computes the instantaneous angular frequency of the variable signal.
% dx represents the spacing between neighbouring points in signal.
% if the spacing is constant, the input dx can be a scalar. Else, it is a vector
% of dimension N - 1.
    [rows, cols] = size(signal);
    if rows>1
        N = rows;
        w = zeros(N, 1);
        if length(dx) == 1
            dx = ones(N-1, 1)*dx;
        end
    else
        N = cols;
        w = zeros(1, N);
        if length(dx) == 1
            dx = ones(1, N-1)*dx;
        end
    end
    signal = hilbert(real(signal)); % + 1i*hilbert(imag(signal));
%     if isreal(signal)
%         signal = hilbert(signal);
%     end
    phi = angle(signal);
    % Making phi a continuous function:
    phi = phi_cont1(phi);
    w(1) = (phi(2) - phi(1))./dx(1);
    w(2:(N - 1)) = (phi(3:N) - phi(1:(N - 2)))./(dx(1:(N - 2)) + dx(2:(N - 1)));
    w(N) = (phi(N) - phi(N - 1))./dx(N - 1);
    w(w>=0) = mod(w(w>=0), 2*pi);
end