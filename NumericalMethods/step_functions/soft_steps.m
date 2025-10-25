function f_smooth = soft_steps(x, fx, xsmooth, soft_fun)
    if nargin<4
        soft_fun = @(x, a, b) softrectfun(x, a, b, 1);
    end
    dim = size(fx, 1);
    f_smooth = zeros(dim, length(xsmooth));
    Nsteps = length(x) - 1;
    for stepi = 1:Nsteps
        f_smooth = f_smooth + fx(:, stepi)*soft_fun(xsmooth, x(stepi), x(stepi + 1));
    end
end