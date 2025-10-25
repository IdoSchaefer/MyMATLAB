function [Jforb, exvals] = fJforb_eq(psi, T, penalforbv)
    Nt = size(psi, 2) - 1;
    dt = T/Nt;
    integw = ones(1, Nt + 1)*dt;
    integw([1, Nt + 1]) = integw([1, Nt + 1])/2;
    exvals = real(sum(conj(psi).*(penalforbv.*psi)));
    Jforb = -sum(integw.*exvals);
end