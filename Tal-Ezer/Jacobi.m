function [x, xtol, Niter] = Jacobi(A, b, x0, tol)
    maxiter = 5000;
    invD = diag(1./diag(A));
    C = eye(size(A, 1)) - invD*A;
    x = x0;
    r = b - A*x;
    iiter = 0;
    while (norm(r)/norm(b)>tol && iiter<=maxiter)
        x = invD*b + C*x;
        r = b - A*x;
        iiter = iiter + 1;
    end
    xtol = norm(r)/norm(b);
    Niter = iiter;
end   