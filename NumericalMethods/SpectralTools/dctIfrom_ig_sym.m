function dctv = dctIfrom_ig_sym(v, interval, int_grid, weights)
% The inverse of the function dctIintgrid; applies for a symmetric internal
% grid structure.
% The function computes the DCT-I of a time signal v, which is
% sampled in a periodic grid structure with a symmetric internal grid in each
% time-step. The time-step length is uniform.
% interval:  The whole interval of v.
% int_grid:  The internal grid inside the interval between equally spaced
% points, in the time space.
% weights: The integration weights of the points in the internal grid.
%%%% dctIfrom_ig_sym1 is more efficient. %%%%
    Nig = length(int_grid);
    dim = size(v);
    if dim(1) == 1
        Np = dim(2);
        v = v.';
    else
        Np = dim(1);
    end
    N = (Np - 1)/Nig + 1;
    Nfourier = 2*(N - 1);
    v_ext = [v(1:Np); v((Np - 1):-1:2)];
    Next = 2*Np - 2;
    ifft_vi = ifft(v_ext(1:Nig:Next));
    dctv = weights(1)*ifft_vi(1:N);
    wdt = (0:2*pi/Nfourier:pi).';
    dt = interval/(N - 1);
    normalig = int_grid/dt;
    for gridi = 2:Nig
        ifft_vi = ifft(v_ext(gridi:Nig:Next));
        dctv = dctv + weights(gridi)*exp(1i*wdt*normalig(gridi)).*ifft_vi(1:N);
    end
    dctv = dctv*sqrt(Nfourier);
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end
end