function W = psi2wigner(psi)
% The function computes the Wigner distribution W from a given wave
% function psi in the x domain.
    N = length(psi);
    % psi should be extended by interpolation, in order to utilize all the
    % data:
    psi_int = interpft(psi, 2*N)/sqrt(2);
    % The density matrix:
    ro = psi_int*psi_int';
    % W is the Wigner distribution matrix:
    W = zeros(2*N);
    % is_in_ro is a logical matrix, which indicates if xi and xtagi are
    % inside the ro matrix. All the x and x' values which are outside ro
    % are set to be 0.
    % ro_v_indices is a matrix which contains the indices of the relevant terms in ro, when is treated as a vector,
    % rather than as a matrix:
    [ro_v_indices, is_in_ro] = get_wigner_indices(N);
    % Here, W, ro an is_in_ro are treated as a vector:
    W(is_in_ro) = ro(ro_v_indices(is_in_ro));
    % Performing a Fourier transform along the s axis. Note that W is
    % shifted, in order that x=0 and s=0 will be in the middle of the axis.
    % It is important to use ifft, with exp(+ipx) basis functions.
    %W = real(ifft(ifftshift(W, 1)));
    W = real(ifft(fftshift(W, 1)));
    % The shift is necessary in order that p=0 will be in the middle of the
    % axis:
    W = fftshift(W, 1);
end