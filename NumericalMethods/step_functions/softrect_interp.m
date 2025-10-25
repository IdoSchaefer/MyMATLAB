function [finterp, coefs] = softrect_interp(x, fxmid, xresult, srect_fun)
    if nargin<4
        srect_fun = @(x, a, b) softrectfun(x, a, b, 1);
    end
    Nsp = length(x) - 1;
    xmid = (x(1:Nsp) + x(2:(Nsp + 1)))/2;
    M = srect_fun(xmid.', x(1:Nsp), x(2:(Nsp + 1)));
    coefs = M\fxmid.';
    finterp = (srect_fun(xresult.', x(1:Nsp), x(2:(Nsp + 1)))*coefs).';
end