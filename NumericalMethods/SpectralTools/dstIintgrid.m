function dstv = dstIintgrid(v, interval, int_grid)
% The function computes the dstI of v (see details in the function dstI), 
% in an internal grid within the small intervals of the equally spaced 
% Fourier grid.
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
    dt = interval/(N + 1);
    Nfourier = 2*(N + 1);
    wdt = (0:2*pi/Nfourier:(2*pi*(1 - 1/Nfourier))).';
    wdt((N + 3):Nfourier) = wdt((N + 3):Nfourier) - 2*pi;
    dstv = zeros(Nfourier*Nig, 1);
    normalig = int_grid/dt;
    for gridi = 1:Nig
        % It is assumed that the function is 0 at the bounaries.
        dstv(gridi:Nig:(Nfourier*Nig)) = 1i/sqrt(Nfourier)*fft([0; exp(-1i*wdt(2:(N + 1))*normalig(gridi)).*v(1:N); ...
             0; -exp(-1i*wdt((N + 3):Nfourier)*normalig(gridi)).*v(N:-1:1)]);
    end    
    dstv = dstv(2:((N + 1)*Nig));
    if dim(1) == 1
        dstv = dstv.';
    end
    if isreal(v)
        dstv = real(dstv);
    elseif isreal(1i*v)
        dstv = 1i*imag(dstv);
    end
end