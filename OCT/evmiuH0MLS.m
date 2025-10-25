function [evmiut, evmiuw, J1, psi] = evmiuH0MLS(psi0, H0, miu, filtermiu, T, dt)
    Nt = T/dt;
    dw = pi/T;
    psi = fMdiag(H0, @(x) exp(-1i*x), psi0, T, Nt);
    evmiut = evmiuMLS(psi, miu);
    dctfactor = T/(sqrt(Nt*pi));
    evmiuw = dctI(evmiut)*dctfactor;
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
    end
    J1fun = 0.5*evmiuw.^2.*vfiltermiu*dw;
    J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
    J1 = sum(J1fun);
end