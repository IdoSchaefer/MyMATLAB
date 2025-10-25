[Vopt649, optval649, flag649, data649] = fminunc(@(V) percosV0lb_grad3(V, [0 64*dx], [0.2 5], 64, 49, 0), Vopt649, options);
Vopt649x = Vopt2Vx0lb3(Vopt649, 64);