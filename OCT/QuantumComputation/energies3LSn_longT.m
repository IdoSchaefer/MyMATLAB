psi0 = [1;0;0];
target = [0;1;0];
H0 = diag([0, 1, 1.96]);
Edomain = [-1, 3];
miu = [0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0];
optionsa = optionsOCqn(1e-4, 1e4);
optionsa.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 0.35);
allTa = (20:10:50)*2*pi*5;
% There's a problem with the convergence of 7*2*pi*5.
NTa = length(allTa);
Nt_tsa = 7;
Ncheba = 7;
Nt_maxa = 1000;
dta = allTa(end)/Nt_maxa;
fieldsa = zeros(NTa, Nt_maxa + 1);
allfieldsa = zeros(NTa, Nt_maxa*(Nt_tsa - 1) + 1);
allpsia = zeros(3, Nt_maxa + 1, length(allTa));
allconva = zeros(NTa, 10001);
allnitera = zeros(1, NTa);
allJ1a = zeros(1, NTa);
energiesa = zeros(1, NTa);
allrelEa = zeros(1, NT);
for Ti = 1:NTa
    Nti = round(Nt_maxa*allTa(Ti)/allTa(end));
    [allfieldsa(Ti, 1:(Nti*(Nt_tsa - 1) + 1)), fieldsa(Ti, 1:(Nti + 1)), allpsia(:, 1:(Nti + 1), Ti), allrelEa(Ti), convi, allnitera(Ti), ~, allJ1a(Ti)] = OCqn(psi0, target, H0, Edomain, miu, @(t) pi/allTa(Ti)*sin(t), 1e-3, optionsa, allTa(Ti), dta, Nt_tsa, Ncheba, 1e-4, 1e4);
    allconva(Ti, 1:(allnitera(Ti) + 1)) = convi;
    energiesa(Ti) = 1e3*(allJ1a(Ti) - convi(allnitera(Ti) + 1));
end