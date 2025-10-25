function Dchi = OCfchi(t, chit, K, Vf, x, field, psi, evmiufil, Nx, Nt, dt)
    tiIpln = RKtiIpln(t, Nt, dt);
    tdata = NewtonIpln4((tiIpln-1)*dt, [psi(:, tiIpln); field(tiIpln); evmiufil(tiIpln)], t);
    psit = tdata(1:Nx);
    field_t = tdata(Nx + 1);
    evmiufil_t = tdata(Nx + 2);
    Dchi = -1i*(ifft(K.*fft(chit)) + (Vf(x) - field_t*x).*chit) - evmiufil_t*x.*psit;
end
