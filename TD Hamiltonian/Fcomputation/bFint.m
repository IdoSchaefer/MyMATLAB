function result = bFint(z, M, order)

    Nz = length(z);
    [sp, wsp] = gaussj(order, 0, 0);
    tic
    integrand = M/2*exp(0.5*z*(1 - sp.')).*(ones(Nz, 1)*(((sp.'+1)./2).^(M-1)));
    result = integrand*wsp;
    toc
end