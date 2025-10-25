function [result, orders] = order_bFint(z, M, tol)
    
    Nz = length(z);
    integrand = zeros(Nz, 100);
    result = zeros(Nz, 1);
    old_result = zeros(Nz, 1);
    orders = zeros(Nz, 1);
    is_not_converged = true(Nz, 1);
    Nz1 = ones(Nz, 1);
%     sp = zeros(100, 1);
%     wsp = zeros(100, 1);
%z=[-100:-1]';  
    order = 2;
    while max(is_not_converged)
        [sp, wsp] = gaussj(order, 0, 0);
        integrand(is_not_converged, 1:order) = M/2*exp(0.5*z(is_not_converged)*(1 - sp.')).*(Nz1(is_not_converged)*(((sp.'+1)./2).^(M-1)));
        old_result(is_not_converged) = result(is_not_converged);
        result(is_not_converged) = integrand(is_not_converged, 1:order)*wsp;
        is_not_converged(is_not_converged) = abs(result(is_not_converged) - old_result(is_not_converged))./abs(result(is_not_converged)) > tol;
        order = order + 1;
        orders(is_not_converged) = order;
    end
    %hold on
    %plot(z,result)
end

