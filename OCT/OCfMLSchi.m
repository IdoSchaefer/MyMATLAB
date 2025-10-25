function Dchi = OCfMLSchi(t, chit, H0, miu, field, psi, evmiufil, Npsi, Nt, dt)
    tiIpln = RKtiIpln(t, Nt, dt);
    tdata = NewtonIpln4((tiIpln-1)*dt, [psi(:, tiIpln); field(tiIpln); evmiufil(tiIpln)], t);
    psit = tdata(1:Npsi);
    field_t = tdata(Npsi + 1);
    evmiufil_t = tdata(Npsi + 2);
    Dchi = -1i*(H0 - field_t*miu)*chit - evmiufil_t*miu*psit;
end
