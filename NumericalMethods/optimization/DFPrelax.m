function [xmin, fmin, gradmin, niter, nfevals, ngevals, reldif_f, reldif_x, alpha] = DFPrelax(f, fgrad, x0, tolf, tolx, maxNiter, alpha0)
% The function minimizes a function by the DFP method. The
% parameter (alpha) which determines the amount of propagation in the
% search direction is determined by a relaxation process.
% Input:
% f: A function handle which represents the function to be optimized
% fgrad: A function handle which represents the gradient function
% x0: The guess solution
% tolf: The tolerance of the function value
% tolx: The tolerance of the solution
% maxNiter: The maximal number of iterations
% alpha0: The guess alpha parameter
% Output:
% xmin: The minimized solution
% fmin: The value of f at xmin
% gradmin: The gradient value at xmin
% niter: The number of required iterations
% nfevals: The number of f evaluations
% ngevals: The number of gradient and inverse Hessian evaluations
% reldif_f: The relative difference of fmin from the f value at the
% previous iteration
% reldif_x: The relative difference of xmin from the solution in the previous iteration
% alpha: The final alpha parameter
    if nargin<7
        alpha = 1;
    else
        alpha = alpha0;
    end
    niter = 0;
    xmin =  x0;
    fmin = f(x0);
    gradmin = fgrad(x0);
    dim = length(x0);
    invHess = eye(dim);
    direction = -gradmin;
    reldif_f = 1 + tolf;
    reldif_x = 1 + tolx;
    nfevals = 1;
    ngevals = 1;
    while (reldif_f>tolf || reldif_x >tolx) && max(gradmin~=0) && niter<maxNiter && alpha>10*eps
        xtrial = xmin + alpha*direction;
        ftrial = f(xtrial);
        nfevals = nfevals + 1;
        % Determining an acceptable alpha value and the corresponding x and
        % f values:
        while ftrial>fmin && alpha>10*eps
            alpha = alpha/2;
            xtrial = xmin + alpha*direction;
            ftrial = f(xtrial);
            nfevals = nfevals + 1;
        end
        reldif_f = abs(ftrial - fmin)/abs(ftrial);
        reldif_x = norm(xtrial - xmin)/norm(xtrial);
        % Making sure that the condition of reldif_f is satisfied due to the fact
        % that we are sufficiently close to the solution; the possibility
        % of a too large alpha should be eliminated.
        if reldif_f<=tolf && alpha>10*eps
            alpha = alpha/2;
            xtrial = xmin + alpha*direction;
            ftrial = f(xtrial);
            nfevals = nfevals + 1;
            reldif_f = abs(ftrial - fmin)/abs(ftrial);
            reldif_x = norm(xtrial - xmin)/norm(xtrial);
        end
        xmin_new = xtrial;
        dx = xmin_new - xmin;
        xmin = xmin_new;
        fmin = ftrial;
        gradmin_new = fgrad(xtrial);
        dgrad = gradmin_new - gradmin;
        gradmin = gradmin_new;
        invHess_dgrad = invHess*dgrad;
        invHess = invHess + (dx*dx.')/(dx.'*dgrad) - (invHess_dgrad*invHess_dgrad.')/(dgrad.'*invHess_dgrad);
        direction = -invHess*gradmin;
        ngevals = ngevals + 1;
        niter = niter + 1;
    end
end