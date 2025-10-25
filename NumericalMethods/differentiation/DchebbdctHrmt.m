function Dv = DchebbdctHrmt(v, domainL)
% The function returns the derivative of a function, multiplied by (Jacobian)^(1/2), sampled at the
% Chebyshev points that include the boundaries of the domain.
% v: A vector that contains the function values at the Chebyshev points.
% domainL: The length of the x domain.
    dim = size(v);
    if dim(1) == 1
        N = dim(2) - 1;
        v = v.';
    else
        N = dim(1) - 1;
    end
    Dv = zeros(N + 1, 1);
    theta = (0:pi/N:pi).';
    k = (0:N).';
%    dctv = dctI([v(1); v(2:N)./(1i*sqrt(sin(theta(2:N)))); v(N+1)]);
    dctv = dctI([0; v(2:N)./sqrt(sin(theta(2:N))); 0]);
    Dv(2:N) = dstI(k(2:N).*dctv(2:N))./sqrt(sin(theta(2:N)))*2/domainL;
%     c = sqrt(2/N)*dctv;
%     c([1, N+1]) = 0.5*c([1, N+1]);
%     sqk = k.^2;
%     Dv(1) = sqk.'*c;
%     sqk_sgn = sqk;
%     sqk_sgn(3:2:(N + 1)) = -sqk_sgn(3:2:(N + 1));
%     Dv(N + 1) = sqk_sgn.'*c;
    if dim(1) == 1
        Dv = Dv.';
    end
end