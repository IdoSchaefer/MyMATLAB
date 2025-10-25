eV = 0.03675;
 beta = 1/(3.1668e-6*298);
Tv = getTv([10 -10 10 -10]*eV, [0 4 12 4], 0, 1, 1e5);
figure
plot((0:1e-5:1)/eV, Tv)
DeltaF1 = sinh(beta*1.29*eV/2)./(cosh(beta*(((0:0.01:1)-5*eV))) + cosh(beta*1.29*eV/2));
hold on
plot((0:1e-2:1)/eV, DeltaF1, 'r')
DeltaF2 = sinh(beta*4.87*eV/2)./(cosh(beta*(((0:0.01:1)-5*eV))) + cosh(beta*4.87*eV/2));
plot((0:1e-2:1)/eV, DeltaF2, 'g')
xlabel('E [eV]')