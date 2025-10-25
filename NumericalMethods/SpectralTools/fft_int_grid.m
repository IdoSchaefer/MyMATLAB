function ftv = fft_int_grid(v, interval, int_grid)
% The function computes the fft of v, in an internal grid within the small
% intervals of the equally spaced Fourier grid.
% interval:  The whole interval of v in the new space (time space).
% int_grid:  The internal grid inside the interval between equally spaced points, in the time space.
    Nig = length(int_grid);
    szv = size(v);
    if szv(1) == 1
        N = szv(2);
        ftv = zeros(1, N*Nig);
        wdt = 0:2*pi/N:(2*pi*(1 - 1/N));
    else
        N = szv(1);
        ftv = zeros(N*Nig, 1);
        wdt = (0:2*pi/N:(2*pi*(1 - 1/N))).';
    end
    dt = interval/N;
    wdt(ceil((N+1)/2):N) = wdt(ceil((N+1)/2):N) - 2*pi;
    normalig = int_grid/dt;
    for gridi = 1:Nig
        ftv(gridi:Nig:(N*Nig)) = fft(exp(-1i*wdt*normalig(gridi)).*v);
    end
end