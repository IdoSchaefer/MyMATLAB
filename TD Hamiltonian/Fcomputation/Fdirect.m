function result = Fdirect(z, t, Nt_ts)
    minus_izt = -1i*z*t;
    Nx = length(z);
    result = exp(minus_izt);
    for polyi = 1:Nt_ts
        result = polyi*(result - 1)./minus_izt;
    end
    result = result.*((ones(Nx, 1)*t).^Nt_ts);
end