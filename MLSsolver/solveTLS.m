% The program assumes the existence of the hamiltonian matrix: H, the final
% time tf,
c0 = [1; 0];
dt = 0.1;
sol = ode45(@TLS, [0 tf], c0, [], H);
c = deval(sol, 0:dt:tf);
P = c.*conj(c);
%figure
plot(0:dt:tf, P(:, :))