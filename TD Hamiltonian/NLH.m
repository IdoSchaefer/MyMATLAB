function Du = NLH(t, u, H0, Vt, x)
    Du = -1i*(H0*u + Vt(u, x, t).*u);
end