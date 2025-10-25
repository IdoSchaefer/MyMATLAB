function [performance, gradient] = percosV_grad(Vk, xdomain, kdomain, Nk, szpenal, pospenal, t)
% The function computes the performance parameter for a cosine series
% absorbing potential. A restriction on the size of the potential can be
% imposed, using a penalty function. The imaginary part can be restricted
% to non-positive values, using an exponential penalty function.
% Vk: A vector that contains the dctI coefficients of the potential. The
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
        Nx = sz(2)/2;
        Vk = Vk.';
    else
        Nx = sz(1)/2;
    end
    Vxreal = real(dctI(Vk(1:Nx)));
    Vximag = real(dctI(Vk((Nx + 1):2*Nx)));
    Vxreim = [Vxreal; Vximag];
    [parameter, grad_norms] = Vabs_efficiency1(Vxreal + 1i*Vximag, xdomain, kdomain, Nk);
    performance = parameter + szpenal*sum(Vxreim.^2) + pospenal*sum(exp(Vximag*t));
    % The gradient of the performance wrt the x potential values (real and
    % imaginary parts separately):
    grad_perf_Vx = grad_norms + 2*szpenal*Vxreim;
    grad_perf_Vx((Nx + 1):2*Nx) = grad_perf_Vx((Nx + 1):2*Nx) + pospenal*t*exp(Vximag*t);
    % The gradient wrt the dctI coefficients:
    gradient = real([dctI([2*grad_perf_Vx(1); grad_perf_Vx(2:(Nx - 1)); 2*grad_perf_Vx(Nx)]);
                     dctI([2*grad_perf_Vx(Nx + 1); grad_perf_Vx((Nx + 2):(2*Nx - 1)); 2*grad_perf_Vx(2*Nx)])]);
    gradient([1; Nx; Nx + 1; 2*Nx]) = gradient([1; Nx; Nx + 1; 2*Nx])/2;  
end