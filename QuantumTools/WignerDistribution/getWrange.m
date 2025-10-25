function Wrange = getWrange(U)
% The function computes the Wigner distribution range of values for many
% time points, for the use of viewWigner.
    szU = size(U);
    Nx = szU(1);
    Nt = szU(2);
    % U should be extended by interpolation, in order to utilize all the
    % data:
    U = interpft(U, 2*Nx)/sqrt(2);
    % The density matrix:
    ro = U(:, 1)*U(:, 1)';
    % W is the Wigner distribution matrix:
    W = zeros(2*Nx);
    % is_in_ro is a logical matrix, which indicates if xi and xtagi are
    % inside the ro matrix. All the x and x' values which are outside ro
    % are set to be 0.
    % ro_v_indices is a matrix which contains the indices of the relevant terms in ro, when is treated as a vector,
    % rather than as a matrix:
    [ro_v_indices, is_in_ro] = get_wigner_indices(Nx);
    % Here, W, ro an is_in_ro are treated as a vector:
    W(is_in_ro) = ro(ro_v_indices(is_in_ro));
    % Performing a Fourier transform along the s axis. Note that W is
    % shifted, in order that x=0 and s=0 will be in the middle of the axis.
    % It is important to use ifft, with exp(+ipx) basis functions.
    W = real(ifft(fftshift(W, 1)));
    % The shift is necessary in order that p=0 will be in the middle of the
    % axis:
    W = fftshift(W, 1);
    Wrange = [min(min(W)), max(max(W))];
    for ti = 2:Nt
        ro = U(:, ti)*U(:, ti)';
        W(is_in_ro) = ro(ro_v_indices(is_in_ro));
%         W = [W((Nx + 1):(2*Nx), :); W(1:Nx, :)];
%         W = real(ifft(W));
        %W = real(ifft([W((Nx + 1):(2*Nx), :); W(1:Nx, :)]));
        W = real(ifft(fftshift(W, 1)));
        W = fftshift(W, 1);
        Wtrange = [min(min(W)), max(max(W))];
        if Wtrange(1)<Wrange(1)
            Wrange(1) = Wtrange(1);
        end
        if Wtrange(2)>Wrange(2)
            Wrange(2) = Wtrange(2);
        end
    end    
end