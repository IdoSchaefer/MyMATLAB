function [evmiuallnt, evmiuallnw, J1alln] = evmiuns(psiE, miu, T, filtermiu, Ndiags)
% The function computes the expectation values of the Hermitian components
% of a Hermitian operator; every component is the n diagonal + (-n)
% diagonal of the matrix.
% psiE: psi in the energy representation, in all time points
% miu: the Hermitian operator
% Ndiags: the number of Hermitian components, that their expectation values at all
% times will be returned.
    if nargin<5
        Ndiags = size(miu, 1);
    end
    Nt = size(psiE, 2) - 1;
    dw = pi/T;
    dctfactor = T/sqrt(Nt*pi);
    vfiltermiu = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfiltermiu(wi) = filtermiu((wi-1)*dw);
    end
    evmiuallnt = zeros(Ndiags, Nt + 1);
    evmiuallnw = zeros(Ndiags, Nt + 1);
    J1alln = zeros(Ndiags, 1);
    evmiuallnt(1, :) = evmiu(psiE, diag(miu));
    evmiuallnw(1, :) = dctI(evmiuallnt(1, :))*dctfactor;
    J1fun = 0.5*evmiuallnw(1, :).^2.*vfiltermiu*dw;
    J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
    J1alln(1) = sum(J1fun);
    for diagi = 1:(Ndiags - 1)
        miui = diag(diag(miu, diagi), diagi) + diag(diag(miu, -diagi), -diagi);
        evmiuallnt(diagi + 1, :) = evmiuMLS(psiE, miui);
        evmiuallnw(diagi + 1, :) = dctI(evmiuallnt(diagi + 1, :))*dctfactor;
        J1fun = 0.5*evmiuallnw(diagi + 1, :).^2.*vfiltermiu*dw;
        J1fun([1, Nt + 1]) = J1fun([1, Nt + 1])/2;
        J1alln(diagi + 1) = sum(J1fun);
    end
end    