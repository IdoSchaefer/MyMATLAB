function Vt = inv_dctIMtp(Vw, t, T, N)
% The function computes the transformation from the frequency domain to the
% time domain of a set of vectors containing the cosine coefficients of
% several functions. The cosine coefficients are represented by the dctI
% representation. The evaluation in the time domain is in arbitrary time points 
% represnted by the vector t.
% Input:
% Vw: A matrix containing the cosine coefficients, where different
% functions are represented by separate columns. It may contain less terms
% than the number of the dctI grid points, if the values of the rest of the
% frequencies is 0.
% t: A vector containing the time-points for evaluation in the time domain.
% Can be either a column vector or a row vector.
% T: The final time
% N: The number of points in the frequency grid minus 1 (the 0 frequency is
% not counted)
% Output:
% Vt: A matrix containing the transformed vectors to the time grid in
% separate columns
    size_t = size(t);
    if size_t(1) == 1
        % If t is a row vector, it is transposed to a column vector:
        t = t.';
    end
    sizeV = size(Vw);
    Nw = sizeV(1) - 1;
    % The omega grid:
    dw = pi/T;
    w = 0:dw:Nw*dw;
    % Computing the transformation matrix:
    transM = cos(t*w);
    transM(:, 1) = transM(:, 1)/2;
    if Nw == N
        transM(:, Nw + 1) = transM(:, Nw + 1)/2;
    end
    Vt = transM*Vw*sqrt(2/N);
end