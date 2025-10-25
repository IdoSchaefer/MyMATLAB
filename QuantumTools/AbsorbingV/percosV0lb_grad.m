function [performance, gradient] = percosV0lb_grad(Vk, xdomain, kdomain, Nk, szpenal, pospenal, t)
% The function computes the performance parameter for a cosine series
% absorbing potential. A zero boundary condition is imposed in the left 
% boundary. A restriction on the size of the potential can be
% imposed, using a penalty function. The imaginary part can be restricted
% to non-positive values, using an exponential penalty function.
% Vk: A vector that contains the dctI coefficients of the potential. 
% The N coefficient is not included. The
% coefficients of the real part are the first values, and of the imaginary
% part are the last. Must be real. 
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nk: The number of equally spaced points in the k domain, (including the
% boundaries).
% szpenal: The penalty parameter on the size of the potential.
% pospenal: A penalty parameter on positive values of the imaginary part.
% It is the coefficient of the exponent.
% t: The coefficient of the exponent argument (time).
    sz = size(Vk);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        Vk = Vk.';
    else
        Nk_free = sz(1)/2;
    end
    Nx = Nk_free + 1;
    % The Vk_N term of the real and imaginary parts:
    VkNre = -2*(Vk(1)/2 + sum(Vk(2:Nk_free)));
    VkNim = -2*(Vk(Nk_free + 1)/2 + sum(Vk((Nk_free + 2):2*Nk_free)));
    Vxreal = real(dctI([Vk(1:Nk_free); VkNre]));
    Vximag = real(dctI([Vk((Nk_free + 1):2*Nk_free); VkNim]));
    Vxreim = [Vxreal; Vximag];
    [parameter, grad_norms] = Vabs_efficiency1(Vxreal + 1i*Vximag, xdomain, kdomain, Nk);
    performance = parameter + szpenal*sum(Vxreim.^2) + pospenal*sum(exp(Vximag(2:Nx)*t));
    % The gradient of the performance wrt the x potential values (real and
    % imaginary parts separately):
    grad_perf_Vx = grad_norms + 2*szpenal*Vxreim;
    grad_perf_Vx((Nx + 2):2*Nx) = grad_perf_Vx((Nx + 2):2*Nx) + pospenal*t*exp(Vximag(2:Nx)*t);
    dctgradVx = real([dctI([2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)]);
                     dctI([2*grad_perf_Vx(Nx + 1); grad_perf_Vx((Nx + 2):(2*Nx - 1)); 2*grad_perf_Vx(2*Nx)])]);
    % The gradient wrt the dctI coefficients:
    gradient = [dctgradVx(1:(Nx - 1)) - dctgradVx(Nx); dctgradVx((Nx + 1):(2*Nx - 1)) - dctgradVx(2*Nx)];
    gradient([1; Nk_free + 1]) = gradient([1; Nk_free + 1])/2;  
end