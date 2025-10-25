load optV40
whos
Vf = @(x)1-1./sqrt(x.^2+1)
xdomain = [-80 80];
xdlenght=160;
Nx=256;
dx=xdlength/Nx
xdlength=160;
clear xdlenght
dx=xdlength/Nx
whos
[Vabs, xabs, x, p, K, Nx] = get_prop_vars(Vf, xdomain, dx, Vopt649x);
[fi0, E0, x, E, P, H] = gsV(Vabs, xdomain, Nx);