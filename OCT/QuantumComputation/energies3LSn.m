psi0 = [1;0;0];
target = [0;1;0];
H0 = diag([0, 1, 1.96]);
Edomain = [-1, 3];
miu = [0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0];
options = optionsOCqn(1e-4, 1e4);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 0.35);
allT = (4:10)*2*pi*5;
% There's a problem with the convergence of 7*2*pi*5.
NT = length(allT);
Nt_ts = 7;
Ncheb = 7;
Nt_max = 1000;
dt = allT(end)/Nt_max;
fields = zeros(NT, Nt_max + 1);
allfields = zeros(NT, Nt_max*(Nt_ts - 1) + 1);
allpsi = zeros(3, Nt_max + 1, length(allT));
allconv = zeros(NT, 10001);
allniter = zeros(1, NT);
allJ1 = zeros(1, NT);
energies = zeros(1, NT);
for Ti = 1:NT
    Nti = round(Nt_max*allT(Ti)/allT(end));
    [allfields(Ti, 1:(Nti*(Nt_ts - 1) + 1)), fields(Ti, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), Ti), ~, convi, allniter(Ti), ~, allJ1(Ti)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.01*sin(t), 1e-3, options, allT(Ti), dt, 7, 7, 1e-4, 1e4);
    allconv(Ti, 1:(allniter(Ti) + 1)) = convi;
    energies(Ti) = 1e3*(allJ1(Ti) - convi(allniter(Ti) + 1));
end