xmax = 200; Nx = 2e4; Emax_eV = 20; NE = 100;
dx = xmax/Nx;
dE_eV = Emax_eV/NE;
x = 0:dx:xmax;
V = 8*0.03675*exp(-x.^2/(2*40^2)).*sin(pi*x/4);
DV = V(2:(Nx + 1)) - V(1:Nx);
Da = ones(1, Nx)*dx;
a(1) = 0;
[Tv checkv]=getTv(DV, Da, 0, Emax_eV*0.03675, NE);
figure
plot(0:dE_eV:Emax_eV, Tv)
title('Periodic barrier');
xlabel('E [eV]')
ylabel('T(E)')