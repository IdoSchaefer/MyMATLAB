function w = instwcos(signal, L)
% The function computes the instantaneous angular frequency of the variable signal.
% Intended to signal which satisfy the cosine series boudary conditions.
% The boundary including grid is used.
% L: The length of the domain.
    signal = hilbert(real(signal)); %+ 1i*hilbert(imag(signal));
%     if isreal(signal)
%         signal = hilbert(signal);
%     end
    phi = angle(signal);
    % Making phi a continuous function:
    phi = phi_cont1(phi);
    w = Dcosines(phi, L);
    w(w>=0) = mod(w(w>=0), 2*pi);
end