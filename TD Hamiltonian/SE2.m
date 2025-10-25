function Du = SE2(t, u, K, Vt, x)
    Du = -1i*(ifft(K.*fft(u)) + Vt(x, t).*u);
end