function [performance, gradient] = percosV0lb_grad7(opt_params, xdomain, kdomain, Nx, Nk, szpenal)
%%%%%%%%%  Applies only for Nk<Nx/2. It should be fixed. %%%%%%%%%%%%
% The function computes the performance parameter for a cosine series
% absorbing potential. A zero boundary condition is imposed in the left 
% boundary. A restriction on the size of the potential can be
% imposed, using a penalty function. The imaginary part is restricted
% to non-positive values, by using an appropriate restriction on the first
% dctI coefficient, and defining the other cefficients as the square of the
% optimized parameters for the imaginary part.
% Vk: A vector that contains the dctI coefficients of the potential. The N
% of the dctI coefficients is Nx, and not the number of elements in Vk.
% The 0 coefficient is not included. The
% coefficients of the real part are the former values, and of the imaginary
% part are the last. Must be real. 
% xdomain: The x domain, [xmin, xmax].
% kdomain: The k domain, [kmin, kmax].
% Nx: The number of x points in the grid minus 1.
% Nk: The number of equally spaced points in the k domain, (including the
% boundaries).
% szpenal: The penalty parameter on the size of the potential.
    sz = size(opt_params);
    if sz(1) == 1
        Nk_free = sz(2)/2;
        opt_params = opt_params.';
    else
        Nk_free = sz(1)/2;
    end
    Vk = [opt_params(1:Nk_free); zeros(2*Nk_free, 1)];
    Vk(Nk_free + (2:2:2*Nk_free)) = opt_params((Nk_free + 1):2*Nk_free).^2;
    Vxreal = dctI0lb_2([Vk(1:Nk_free); zeros(Nx - Nk_free, 1)]);
    Vximag = [0; dstI([Vk((Nk_free + 1):3*Nk_free); zeros(Nx - 2*Nk_free - 1, 1)]); 0] -...
        sqrt(2/Nx^3)*pi*sum(((2:2:2*Nk_free).').*Vk((Nk_free + 2):2:3*Nk_free))*((0:Nx).'); 
    Vxreim = [Vxreal; Vximag];
    [parameter, grad_norms] = Vabs_efficiency1(Vxreal + 1i*Vximag, xdomain, kdomain, Nk);
    performance = parameter + szpenal*sum(Vxreim.^2);
    % The gradient of the performance wrt the x potential values (real and
    % imaginary parts separately):
    grad_perf_Vx = grad_norms + 2*szpenal*Vxreim;
    transgradVx = real([dctI([2*grad_perf_Vx(1); grad_perf_Vx(2:Nx); 2*grad_perf_Vx(Nx + 1)]);...
                     dstI(grad_perf_Vx((Nx + 3):(2*Nx + 1)))]);
    % The gradient wrt the dctI coefficients:
    gradient = [transgradVx(2:(Nk_free + 1)) - transgradVx(1); ...
        2*opt_params((Nk_free + 1):2*Nk_free).*(transgradVx((Nx + 3):2:(Nx + 1 + 2*Nk_free)) -...
                     sqrt(2/Nx^3)*pi*sum(((1:Nx).').*grad_perf_Vx((Nx + 3):2*(Nx + 1)))*((2:2:2*Nk_free).'))];
end