function [params, Vk] = Vx2params0lb5(Vx, Nk_free)
% The function converts the x potential values for a 0 left boundary to the
% dctI coefficients, for the use of percosV0lb_grad5.
% Vx includes the V(x) values, without the 0 boundary, for the real and
% imaginary parts separately.
    dim = size(Vx);
    if dim(1) == 1
        N = dim(2)/2;
        Vx = Vx.';
    else
        N = dim(1)/2;
    end
    Vk = [dctI([0; Vx(1:N)]); dctI([0; Vx((N + 1):2*N)])];
    % This is just one possibility for the imaginary params:
    params = real([Vk(2:(Nk_free + 1)); sqrt(Vk((N + 3):(N + Nk_free + 2)))]);
    if dim(1) == 1
        params = params.';
        Vk = Vk.';
    end    
end