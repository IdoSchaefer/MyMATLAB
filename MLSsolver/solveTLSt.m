% The program assumes the existence of the vector of the hamiltonian eigenvalues: E,
% the functional form of the time dependent purturbation, @coupling(t), the final time tf,
c0 = [1; 0];
dt = 0.1;
sol = ode45(@TLSt, [0 tf], c0, [], E, coupling);
c = deval(sol, 0:dt:tf);
P = c.*conj(c);
%figure
plot(0:dt:tf, P(:, :))