function U = fMdiag(M, f, u0, T, Nt)
% The program computes: u(t) = f(M*t)*u0, when M is a matrix, and u0 is a vector. The
% result is computed at Nt equally spaced time points, when T is the maximal
% time. U contains the vector results, in different columns for different
% time points.
% f is a function handle.
    dt = T/Nt;
    dim = length(u0);
    U = zeros(dim, Nt + 1);
    U(:, 1) = u0;
    [P, D] = eig(M);
    u0Dbase = P\u0;
    eigval = diag(D);
    for ti = 1:Nt
        U(:, ti + 1) = f(eigval*ti*dt).*u0Dbase;
    end
    U(:, 2:Nt+1) = P*U(:, 2:Nt+1);
end