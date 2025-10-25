function Dv = DfourierMap(v, r) %domainL)
% The function returns the derivative of a periodic function sampled at the 
% points r, using mapped Fourier differentiation.
% In the current version, there is no advantage over a simple finite
% difference method.
% v: A vector that contains the function values. It should satisfy periodic boundary conditions. 
% r: A vector that contains the sampling points.
%%%%%% domainL: The length of the domain.
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
        r = r.';
    else
        N = dim(1);
    end
%     dx = domainL/N;
%     k = (0:(2*pi/domainL):(2*pi*(1/dx - 1/domainL))).';
%     k((N/2 + 1):N) = k((N/2 + 1):N) - 2*pi/dx;
% %     J = Dcosines(r, domainL);
%     J = zeros(N, 1);
%     J(1) = (r(2) - r(1))/dx;
%     J(N) = (r(N) - r(N - 1))/dx;
%     J(2:(N - 1)) = (r(3:N) - r(1:(N - 2)))./(2*dx);
    k = (0:(2*pi/N):(2*pi*(N - 1)/N)).';
    k((N/2 + 1):N) = k((N/2 + 1):N) - 2*pi;
    J = zeros(N, 1);
    J(1) = (r(2) - r(1));
    J(N) = (r(N) - r(N - 1));
    J(2:(N - 1)) = (r(3:N) - r(1:(N - 2)))/2;
    Dv = ifft(1i*k.*fft(v))./J(1:N);
%    Dv = Dcosines(v, domainL)./Dcosines(r, domainL);
    if dim(1) == 1
        Dv = Dv.';
    end
    if dim(1) == 1
        Dv = Dv.';
    end
end