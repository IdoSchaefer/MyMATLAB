function [performance, gradient] = percosV0lb_grad2(Vk, xdomain, kdomain, Nxpenal, Nk, szpenal, pospenal, t)
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
    %Nx = Nk_free + 1;
    % The Vk_N term of the real and imaginary parts:
    VkNre = -2*(Vk(1)/2 + sum(Vk(2:Nk_free)));
    VkNim = -2*(Vk(Nk_free + 1)/2 + sum(Vk((Nk_free + 2):2*Nk_free)));
    Vxreal = real(dctI([Vk(1:Nk_free); VkNre]));
    Vximag = real(dctI([Vk((Nk_free + 1):2*Nk_free); VkNim]));
    Nxfactor = sqrt(Nxpenal/Nk_free);
    if Nxpenal>Nk_free
        Vximag_penal = dctI([Vk((Nk_free + 1):2*Nk_free); VkNim/2; zeros(Nxpenal - Nk_free, 1)]*Nxfactor);
    else
        Vximag_penal = dctI0lb(Vk((Nk_free + 1):2*Nk_free));
    end    
    Vxreim = [Vxreal; Vximag];
    [parameter, grad_norms] = Vabs_efficiency1(Vxreal + 1i*Vximag, xdomain, kdomain, Nk);
    performance = parameter + szpenal*sum(Vxreim.^2) + pospenal*sum(exp(Vximag_penal(2:(Nxpenal + 1))*t));
    % The gradient of the performance wrt the x potential values (real and
    % imaginary parts separately):
    grad_perf_Vx = grad_norms + 2*szpenal*Vxreim;
    grad_pospenal = [0; pospenal*t*exp(Vximag_penal(2:(Nxpenal + 1))*t)];
    dctgradVx = real([dctI([2*grad_perf_Vx(1); grad_perf_Vx(2:Nk_free); 2*grad_perf_Vx(Nk_free + 1)]);
          dctI([2*grad_perf_Vx(Nk_free + 2); grad_perf_Vx((Nk_free + 3):(2*Nk_free + 1)); 2*grad_perf_Vx(2*(Nk_free + 1))])]);
    dctgrad_pospenal = dctI([0; grad_pospenal(2:Nxpenal); 2*grad_pospenal(Nxpenal + 1)]*Nxfactor);
    % The gradient wrt the dctI coefficients:
    gradient = [dctgradVx(1:Nk_free) - dctgradVx(Nk_free + 1); 
        dctgradVx((Nk_free + 2):(2*Nk_free + 1)) - dctgradVx(2*(Nk_free + 1)) + ...
        dctgrad_pospenal(1:Nk_free) - dctgrad_pospenal(Nk_free + 1)];
    gradient([1; Nk_free + 1]) = gradient([1; Nk_free + 1])/2;  
end