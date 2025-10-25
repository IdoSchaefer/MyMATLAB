function Vk = Vx2Vk0lb(Vx)
% The function converts the x potential values for a 0 left boundary to the
% dctI coefficients, for the use of percosV0lb_grad.
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
    Vk = Vk([1:N, (N + 2):(2*N + 1)]);
    if dim(1) == 1
        Vk = Vk.';
    end    
end