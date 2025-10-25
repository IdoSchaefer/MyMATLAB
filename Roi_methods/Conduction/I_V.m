eV = 0.03675;
Temp = 298;
boundaries = [0.01 0.645 0.655 2.53 2.61 5.1 6.6 8.4]*eV;
allh = [1e-2 1e-6 1e-2 1e-5 1e-2 1e-5 1e-3]*eV;
allL = boundaries(2:end) - boundaries(1:(end - 1));
allNpoints = round(allL./(allh*2))*2;
allh_correct = allL./allNpoints;
allE = cell(7, 1);
allTv = cell(7, 1);
hV = 1e-3;
NV = 6/hV;
allV = (0:hV:6)*eV;
allI = zeros(1, NV + 1);
for Li = 1:7
    allE{Li} = boundaries(Li):allh_correct(Li):boundaries(Li + 1);
    allTv{Li} = getTv([10 -10 10 -10]*eV, [0 4 12 4], boundaries(Li), boundaries(Li + 1), allNpoints(Li));
end
beta = 1/(3.1668e-6*Temp);
for Vi = 1:(NV + 1)
    for Li = 1:7
        DeltaF_Li = sinh(beta*allV(Vi)/2)./(cosh(beta*((allE{Li}-5*eV))) + cosh(beta*allV(Vi)/2));
        allI(Vi) = allI(Vi) + eq_space_simpson(DeltaF_Li.*allTv{Li}, boundaries(Li), boundaries(Li + 1))/pi;
    end
end
figure
plot(allV/eV, allI)
xlabel('V [Volts]')
ylabel('I(V) [a.u.]')