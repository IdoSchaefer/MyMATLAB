eV = 0.03675;
[Tv checkv]=getTv([10 -10 10 -10]*eV, [0 4 12 4], 0, 1, 1e5);
figure
plot((0:1e-5:1)/eV, Tv)
title('Double barrier');
xlabel('E [eV]')
ylabel('T(E)')