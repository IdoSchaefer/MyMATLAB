% The program assumes the existence of the vector of the hamiltonian eigenvalues: E,
% the functional form of the time dependent purturbation, @coupling(t), the final time tf,
c0 = zeros(length(E), 1);
c0(1) = 1;
dt = 0.1;
sol = ode45(@MLS, [0 tf], c0, [], E, coupling);
c = deval(sol, 0:dt:tf);
P = c.*conj(c);
%figure
plot(0:dt:tf, P(:, :))