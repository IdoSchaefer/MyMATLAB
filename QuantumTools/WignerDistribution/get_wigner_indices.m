function [ro_v_indices, is_in_ro] = get_wigner_indices(N)
% The function returns the indices of the relevant terms in the density
% matrix for the use of psi2wigner, and the terms which are inside ro.
% N: The size of the original x domain, before interpolating.
    % The indices of the q and s coordinates in W (before Fourier
    % transformig s):
    qi = ones(2*N, 1)*(1:2*N);
    si = (1:2*N).'*ones(1, 2*N);
    % The indices of x and x' in ro, in the terms of the q = (x + x')/2 and
    % s = x - x' indices in W:
%     xi = qi + si - N;
%     xtagi = qi - si + N;
% si is replaced by si-1, in order to shift all the rows a single row
% downward. This makes the DFT symmetric arround the (2*N)/2 + 1 point.
    xi = qi + si - 1 - N;
    xtagi = qi - si + 1 + N;
    % is_in_ro is a logical matrix, which indicates if xi and xtagi are
    % inside the ro matrix. All the x and x' values which are outside ro
    % are set to be 0.
    is_in_ro = xi>=1 & xi<=2*N & xtagi>=1 & xtagi<=2*N;
    % ro_v_indices is a matrix which contains the indices of the relevant terms in ro, when is treated as a vector,
    % rather than as a matrix:
    ro_v_indices = (xi - 1)*(2*N) + xtagi;
end