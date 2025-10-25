function [performance, gradient] = perVabs(Vxreim, xdomain, kdomain, Nk, szpenal)

    Nx = length(Vxreim)/2;
    [parameter, grad_norms] = Vabs_efficiency1(Vxreim(1:Nx) + 1i*Vxreim((Nx + 1):2*Nx), xdomain, kdomain, Nk);
    performance = parameter + szpenal*sum(Vxreim.^2);
    gradient = grad_norms + 2*szpenal*Vxreim;
end