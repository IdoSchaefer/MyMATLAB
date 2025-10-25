[Tv checkv]=getTv([0.3675 -0.3675], [0 5/0.52918], 0, 1, 100);
figure
plot((0:0.01:1)/0.03675, Tv)
title('Rectangular barrier');
xlabel('E [eV]')
ylabel('T(E)')