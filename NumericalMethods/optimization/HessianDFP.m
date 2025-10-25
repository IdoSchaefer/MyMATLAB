function invHess_new = HessianDFP(invHess, dx, dgrad)
% The function computes the updated inverse Hessian in the DFP method.
% invHess: The previous inverse Hessian
% dx: The difference of the new solution vector from the old one (delta in Fletcher)
% dgrad: The difference of the new gradient vector from the old one (g in Fletcher)
    invHess_dgrad = invHess*dgrad;
    invHess_new = invHess + (dx*dx.')/(dx.'*dgrad) - (invHess_dgrad*invHess_dgrad.')/(dgrad.'*invHess_dgrad);
end