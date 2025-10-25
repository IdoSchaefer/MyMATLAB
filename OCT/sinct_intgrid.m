function sinctv = sinct_intgrid(v, interval, int_grid)
% The function computes a discrete sinc transform of a vector v in an 
% internal grid within the small intervals of the equally spaced Fourier grid.
% The transform without the internal grid is defined in the following way:
% u_k = sqrt(2/pi)*sum_{j=0}^{N - 1}(1/d_j)*v_j*sin(j*k*pi/(N - 1))/j
% interval:  The whole interval of v in the new space (time space).
% int_grid:  The internal grid inside the interval between equally spaced points.
    Nig = length(int_grid);
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
    else
        N = dim(1);
    end
    if size(int_grid, 1) == 1
        int_grid = int_grid.';
    end
    dt = interval/(N - 1);
    dw = pi/interval;
    Nfourier = 2*(N - 1);
    wdt = (0:2*pi/Nfourier:(2*pi*(1 - 1/Nfourier))).';
    wdt((N + 1):Nfourier) = wdt((N + 1):Nfourier) - 2*pi;
    sinctv = zeros(Nfourier*Nig, 1);
    normalig = int_grid/dt;
    for gridi = 1:Nig
        sinctv(gridi:Nig:(Nfourier*Nig)) = 1/sqrt(2*pi)*fft([0; 1i*exp(-1i*wdt(2:(N - 1))*normalig(gridi)).*v(2:(N - 1))./(1:(N - 2)).'; ...
             sin(pi*normalig(gridi))*v(N)/(N - 1); -1i*exp(-1i*wdt((N + 1):Nfourier)*normalig(gridi)).*v((N - 1):-1:2)./((N - 2):-1:1).']);
         % The w=0 term will be added afterwards.
    end    
    sinctv = sinctv(1:((N - 1)*Nig + 1));
    % Adding the w=0 term:
    sinctv = sinctv + 1/sqrt(2*pi)*[kron(wdt(1:(N - 1)), ones(Nig, 1)) + kron(ones((N - 1), 1)*wdt(2), normalig); pi];
    if dim(1) == 1
        sinctv = sinctv.';
    end
    if isreal(v)
        sinctv = real(sinctv);
    elseif isreal(1i*v)
        sinctv = 1i*imag(sinctv);
    end
end