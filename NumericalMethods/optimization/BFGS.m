function [xmin, fmin, gradmin, niter, nfevals, reldif_f, reldif_x] = BFGS(f, fgrad, x0, invHess0, tolf, tolx, maxNiter, sigma, Deltaf0, minimal_f)
% The function minimizes a function by the BFGS method.
% Input:
% f: A function handle which represents the function to be optimized
% fgrad: A function handle which represents the gradient function
% x0: The guess solution
% invHess0: The initial approximated inverse Hessian matrix; when an
% approximation for the Hessian is unavailable, put [] instead to use the identity matrix.
% tolf: The tolerance of the function value
% tolx: The tolerance of the solution
% maxNiter: The maximal number of iterations
% sigma: A parameter which determines the accuracy of the line-search; see
% in Fletcher's book. Put [] instead for the default value (0.9).
% Deltaf0: The estimated reduction in f in the first iteration. If it is
% unknown, put [] instead.
% minimal_f: Defined as \bar{f} in the Fletcher's book; if irrelevat, substitute -inf.
% Output:
% xmin: The minimized solution
% fmin: The value of f at xmin
% gradmin: The gradient value at xmin
% niter: The number of required iterations
% nfevals: The number of f and fgrad evaluations
% reldif_f: The relative difference of fmin from the f value at the
% previous iteration
% reldif_x: The relative difference of xmin from the solution in the previous iteration
    if nargin<8
        sigma = 0.9;
    end
    if nargin<10
        minimal_f = -inf;
    end
    niter = 0;
    xmin_old =  x0;
    fmin_old = f(x0);
    gradmin_old = fgrad(x0);
    if isempty(invHess0)
        dim = length(x0);
        invHess = eye(dim);
        direction = -gradmin_old;
    else
        invHess = invHess0;
        direction = -invHess*gradmin_old;
    end
    
    reldif_f = 1 + tolf;
    reldif_x = 1 + tolx;
    nfevals = 1;
    sol_achieved = false;
    if nargin<9 || isempty(Deltaf0)
        alpha1 = 1;
    else
        alpha1 = min([1, -2*Deltaf0/(gradmin_old.'*direction)]);
    end
    maxNiter_ls = 20;
    while (reldif_f>tolf || reldif_x>tolx) && max(gradmin_old~=0) && niter<maxNiter && ~sol_achieved
%         if niter==56
%             keyboard
%         end
        [xmin, fmin, gradmin, sol_achieved, nfevals_ls] = line_search(f, fgrad, xmin_old, fmin_old, gradmin_old, direction,...
            alpha1, 0.01, sigma, 9, 0.1, 0.5, minimal_f, maxNiter_ls);
        nfevals = nfevals + nfevals_ls;
        dx = xmin - xmin_old;
        dgrad = gradmin - gradmin_old;
        dx_dot_dgrad = dx.'*dgrad;
        invHess_dgrad = invHess*dgrad;
        dx_scaledT = dx.'/dx_dot_dgrad;
        %invHess_dgrad_dxT = invHess_dgrad*dx.';
        invHess_dgrad_dxsT = invHess_dgrad*dx_scaledT;
        %invHess = invHess + ((1 + (dgrad.'*invHess_dgrad)/dx_dot_dgrad)*(dx*dx.') - (invHess_dgrad_dxT + invHess_dgrad_dxT.'))/dx_dot_dgrad;
        invHess = invHess + ((1 + (dgrad.'*invHess_dgrad)/dx_dot_dgrad)*(dx*dx_scaledT) - (invHess_dgrad_dxsT + invHess_dgrad_dxsT.'));
        direction = -invHess*gradmin;
        reldif_f = 2*abs(fmin - fmin_old)/(abs(fmin) + abs(fmin_old));
        reldif_x = 2*norm(xmin - xmin_old)/(norm(xmin) + norm(xmin_old));
        Deltaf = max([fmin_old - fmin, 10*eps]);
        alpha1 = min([1, -2*Deltaf/(gradmin.'*direction)]);
        xmin_old = xmin;
        fmin_old = fmin;
        gradmin_old = gradmin;
        maxNiter_ls = min([3*nfevals_ls, 5]);
        niter = niter + 1;
    end
    if sol_achieved
        fprintf('\nThe optimization process has achieved the minimal function value.\n')
    end
end