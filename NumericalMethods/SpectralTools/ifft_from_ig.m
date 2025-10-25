function iftv = ifft_from_ig(v, interval, int_grid, weights)
% The inverse of the function fft_int_grid;
% The function computes the inverse DFT of a time signal v, which is
% sampled in a periodic grid structure.
% interval:  The whole interval of v, including the upper boundary.
% int_grid:  The internal grid inside the interval between equally spaced
% points, in the time space.
% weights: The integration weights of the points in the internal grid.
    Nig = length(int_grid);
    szv = size(v);
    if szv(1) == 1
        N = szv(2)/Nig;
        iftv = zeros(1, N);
        wdt = 0:2*pi/N:(2*pi*(1 - 1/N));
    else
        N = szv(1)/Nig;
        iftv = zeros(N, 1);
        wdt = (0:2*pi/N:(2*pi*(1 - 1/N))).';
    end
    wdt(ceil((N+1)/2):N) = wdt(ceil((N+1)/2):N) - 2*pi;
    dt = interval/N;
    normalig = int_grid/dt;
    for gridi = 1:Nig
        iftv = iftv + weights(gridi)*exp(1i*wdt*normalig(gridi)).*ifft(v(gridi:Nig:(N*Nig)));
    end
end