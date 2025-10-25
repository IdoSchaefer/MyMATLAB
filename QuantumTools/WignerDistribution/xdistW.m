function [xdist, xint] = xdistW(U, xdomain, pfun)
    szU = size(U);
    Nx = szU(1);
    Nt = szU(2);
    xdlength = xdomain(2) - xdomain(1);
    % U should be extended by interpolation, in order to utilize all the
    % data:
    U = interpft(U, 2*Nx)/sqrt(2);
    dx = xdlength/Nx;
    dp = 2*pi/xdlength;
    maxP = pi/dx;
    xint = xdomain(1):dx/2:(xdomain(2) - dx/2);
    pint = (-maxP:dp/2:(pi/dx - dp/2)).';
%     vx = xfun(xint);
%     maxvx = max(abs(vx));
%     is_none0 = abs(vx)/maxvx>eps;
%     xM = ones(2*Nx, 1)*vx(is_none0);
%     is_none0M = true(2*Nx, 1)*is_none0;
    % The matrix which represents the x function:
    pM = pfun(pint)*ones(1, 2*Nx);
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
    % Here, W, ro and is_in_ro are treated as a vector:
    W(is_in_ro) = ro(ro_v_indices(is_in_ro));
    % Performing a Fourier transform along the s axis. Note that W is
    % shifted, in order that x=0 and s=0 will be in the middle of the axis.
    % It is important to use ifft, with exp(+ipx) basis functions.
    W = real(ifft(fftshift(W, 1)));
    % The shift is necessary in order that p=0 will be in the middle of the
    % axis:
    W = fftshift(W, 1);
    xdist = sum(W.*pM);
    % pdist = sum(W(is_none0M).*xM), 2);
    for ti = 2:Nt
        ro = U(:, ti)*U(:, ti)';
        W(is_in_ro) = ro(ro_v_indices(is_in_ro));
        W = real(ifft(fftshift(W, 1)));
        W = fftshift(W, 1);
        xdist = xdist + sum(W.*pM);    
    end    
end