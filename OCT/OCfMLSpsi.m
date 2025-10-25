function Dpsi = OCfMLSpsi(t, psit, H0, miu, field, Nt, dt)
    tiIpln = RKtiIpln(t, Nt, dt);
    field_t = NewtonIpln4((tiIpln-1)*dt, field(tiIpln), t);
    Dpsi = -1i*(H0 - field_t*miu)*psit;
end
