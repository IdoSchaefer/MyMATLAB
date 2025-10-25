function U = fMkrsimp(M, f, u0, T, Nt, Nkr)
% The program computes: u(t) = f(M*t)*u0, where M is a matrix, and u0 is a vector. The
% result is computed at Nt equally spaced time points, where T is the maximal
% time. U contains the vector results, in different columns for different
% time points.
% The simple Krylov algorithm is used.
% f is a function handle.
    dt = T/Nt;
    dim = length(u0);
    U = zeros(dim, Nt + 1);
    U(:, 1) = u0;
    [V, H] = createKr(M, u0, Nkr);
    [P, D] = eig(H(1:Nkr, :));
    u0kr = zeros(Nkr, 1);
    u0kr(1) = norm(u0);
    u0Dbase = P\u0kr;
    eigval = diag(D);
    U(:, 2:(Nt + 1)) = V(:, 1:Nkr)*P*spdiags(u0Dbase, 0, Nkr, Nkr)*f(eigval*(dt:dt:T));
end 