function [Vabs, xabs, x, p, K, Nx, V0] = get_prop_vars(Vf, xdomain, dx, Vedge)
minx = xdomain(1);
maxx = xdomain(2);
xdlength = maxx - minx;
x = (minx:dx:(maxx-dx)).';
Nx = xdlength/dx;
Nabs = length(Vedge);
edgeL = (Nabs - 1)*dx;
xabs = bound0der1(@(x) x, x, minx + edgeL + 2.5, maxx - edgeL - 2.5, 1);
V0 = bound0der1(Vf, x, minx + edgeL + 2.5, maxx - edgeL - 2.5, 1);
Vabs = [V0(1:Nabs) + Vedge(Nabs:-1:1); V0((Nabs + 1):(Nx - Nabs + 1)); V0((Nx - Nabs + 2):Nx) + Vedge(1:(Nabs - 1))];
p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
end