eV = 0.03675;
Temp = 298;
beta = 1/(3.1668e-6*Temp);
hV = 1e-3;
NV = 6/hV;
hVg = 1e-1;
NVg = 6/hVg;
allV = (0:hV:6)*eV;
allVg_eV = (-3:hVg:3);
allg = zeros(NVg + 1, NV + 1);
hE = 1e-3;
NE = (12 - 0.01)/hE;
E = (0.01:hE:12)*eV;
for Vgi = 1:(NVg + 1)
    T_Vg = getTv([10, -10 + allVg_eV(Vgi), 10 - allVg_eV(Vgi), -10]*eV, [0 4 12 4], 0.01*eV, 12*eV, NE);
    for Vi = 1:(NV + 1)
        VEfun = (cosh(beta*allV(Vi)/2)*cosh(beta*(E - 5*eV)) + 1)./(cosh(beta*(E - 5*eV)) + cosh(beta*allV(Vi)/2)).^2;
        allg(Vgi, Vi) = beta/(2*pi)*eq_space_simpson(VEfun.*T_Vg, 1e-2*eV, 12*eV);
    end
end
figure
surfl(allV/eV,allVg_eV, allg)
shading interp
colormap copper