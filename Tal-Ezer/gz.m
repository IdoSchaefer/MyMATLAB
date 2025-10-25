function [x, xtol, Niter] = gz(A, b, x0, tol)
    maxiter = 5000;
    invL = inv(tril(A));
    C = eye(size(A, 1)) - invL*A;
    x = x0;
    r = b - A*x;
    iiter = 0;
    while (norm(r)/norm(b)>tol && iiter<=maxiter)
        x = invL*b + C*x;
        r = b - A*x;
        iiter = iiter + 1;
    end
    xtol = norm(r)/norm(b);
    Niter = iiter;
end   