function Dpsi = OCguess(t, psi, K, Vf, x, field)
    Dpsi = -1i*(ifft(K.*fft(psi)) + (Vf(x) - field(t)*x).*psi);
end
