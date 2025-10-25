function Du = SE(t, u, H0, Vt, x)
    Du = -1i*(H0*u + Vt(x, t).*u);
end