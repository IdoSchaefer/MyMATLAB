[U mniter matvecs] = TDHmatKr([0 0; 0 0], @(u,t) [0, -1; -1, 0]*0.5e-3*sin(pi*t/9000)^2, [], [1; 0], [0 9e3], 100, 8, 2, 1e-5);
plot(0:90:9e3, conj(U).*U)