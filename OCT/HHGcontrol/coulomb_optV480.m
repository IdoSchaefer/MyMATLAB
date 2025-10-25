[Vabs480, xabs480, x480, p480, K480, Nx480, V0480] = get_prop_vars(Vf, [-480 480], 40/64, Vopt649x);
[fi0480, E0480, ~, E480, ~, ~] = gsV(Vf, [-480 480], 1536);
save coulomb_optV480 Vabs480 xabs480 x480 p480 K480 Nx480 V0480 Vf Vopt649x fi0480 E0480 E480
