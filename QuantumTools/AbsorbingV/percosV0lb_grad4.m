function [performance, gradient] = percosV0lb_grad4(Vk, xdomain, kdomain, Nx, Nk, szpenal)
% The function computes the performance parameter for a cosine series
% absorbing potential. A zero boundary condition is imposed in the left 
% boundary. A restriction on the size of the potential can be
% imposed, using a penalty function. The imaginary part is restricted
% to non-positive values, by defining it as minus the square of the cosine
% series.
% Vk: A vector that contains the dctI coefficients of the potential. The N
% of the dctI coefficients is Nx, and not the number of elements in Vk.
% The N coefficient is not included. The
% coefficients of the real part are the former values, and of the imaginary
% part are the last. Must be real. 
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nx: The number of x points in the grid minus 1.
% Nk: The number of equally spaced points in the k domain, (including the
% boundaries).
% szpenal: The penalty parameter on the size of the potential.
    sz = size(Vk);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        Vk = Vk.';
    else
        Nk_free = sz(1)/2;
    end
    %Nxfactor = sqrt(Nx/Nk_free);
    if Nx>Nk_free        
        Vxreal = dctI([Vk(1:Nk_free); -(Vk(1)/2 + sum(Vk(2:Nk_free))); zeros(Nx - Nk_free, 1)]); %*Nxfactor);
    else
        Vxreal = dctI0lb(Vk(1:Nk_free));
    end    
    dctVkimag = dctI([Vk((Nk_free + 1):2*Nk_free); -(Vk(Nk_free + 1)/2 + sum(Vk((Nk_free + 2):2*Nk_free))); zeros(2*Nx - Nk_free, 1)]); %*Nxfactor); 
    Vximag = -dctVkimag(1:(Nx + 1)).^2;
    Vxreim = [Vxreal; Vximag];
    [parameter, grad_norms] = Vabs_efficiency1(Vxreal + 1i*Vximag, xdomain, kdomain, Nk);
    performance = parameter + szpenal*sum(Vxreim.^2);
    % The gradient of the performance wrt the x potential values (real and
    % imaginary parts separately):
    grad_perf_Vx = grad_norms + 2*szpenal*Vxreim;
    dctgradVx = real([dctI([2*grad_perf_Vx(1); grad_perf_Vx(2:Nx); 2*grad_perf_Vx(Nx + 1)]);
                     dctI([-2*dctVkimag(1:(Nx + 1)).*[2*grad_perf_Vx(Nx + 2); grad_perf_Vx((Nx + 3):(2*(Nx + 1)))]; zeros(Nx ,1)])]);
    % The gradient wrt the dctI coefficients:
    gradient = [dctgradVx(1:Nk_free) - dctgradVx(Nk_free + 1); dctgradVx((Nx + 2):(Nx + 1 + Nk_free)) - dctgradVx(Nx + 2 + Nk_free)];
    gradient([1; Nk_free + 1]) = gradient([1; Nk_free + 1])/2;  
end