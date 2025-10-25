function Dpsi = OCfpsi(t, psit, K, Vf, x, field, Nt, dt)
    tiIpln = RKtiIpln(t, Nt, dt);
    field_t = NewtonIpln4((tiIpln-1)*dt, field(tiIpln), t);
    Dpsi = -1i*(ifft(K.*fft(psit)) + (Vf(x) - field_t*x).*psit);
end
