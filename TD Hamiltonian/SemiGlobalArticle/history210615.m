L = 16*sqrt(pi);
Nx = 128;
dx = L/Nx;
x = (-L/2:dx:(L/2 - dx)).';
% Constructing the kinetic energy matrix diagonal in the p domain:
p = (0:(2*pi/L):(2*pi*(1/dx - 1/L))).';
p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
K = p.^2/2;
% The potential energy matrix diagonal in the x domain:
V = x.^2/2;
v0 = zeros(128,1);
v0(1:10) = 1/sqrt(10);
norm(v0)
[V,H] = createKrop(@(psi) Hpsi(K,V,psi), v0, 10);
eig(H(1:10,:))
ev = eig(H(1:10,:));
ev=sort(ev)
Vmat = diag(V);
% The kinetic energy matrix in the x domain:
Kmat = Nx*ifft(ifft(diag(K))')';
% The Hamiltonian:
H = Kmat + Vmat;
V = x.^2/2;
Vmat = diag(V);
% The kinetic energy matrix in the x domain:
Kmat = Nx*ifft(ifft(diag(K))')';
% The Hamiltonian:
H = Kmat + Vmat;
evfull = sort(eig(H));
evfull
[P, D] = eig(H);
E = diag(D);
%    [E, orderE] = sort(E);
[Ereal, orderE] = sort(real(E));
E = E(orderE);
P = P(:, orderE);
v0 = zeros(128,1);
v0 = 1/sqrt(10)*sum(P(:,1:10),2);
size(v0)
norm(v0)
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 10);
eig(Hess)
eig(Hess(1:10,:))
sort(ans)
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 9);
sort(eig(Hess))
size(Hess)
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 11);
size(Hess)
sort(eig(Hess(1:11,:)))
v0 = 1/sqrt(5)*sum(P(:,1:5),2);
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 5);
sort(eig(Hess(1:5,:)))
figure
plot(x, conj(P(:,5)).*p(:,5))
plot(x, conj(P(:,5)).*P(:,5))
plot(x, conj(P(:,10)).*P(:,10))
x(1)
8*sqrt(pi)
v0 = 1/sqrt(7)*sum(P(:,1:7),2);
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 7);
sort(eig(Hess(1:5,:)))
sort(eig(Hess(1:7,:)))
v0 = 1/sqrt(9)*sum(P(:,1:9),2);
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 9);
sort(eig(Hess(1:7,:)))
sort(eig(Hess(1:9,:)))
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 8);
v0 = 1/sqrt(8)*sum(P(:,1:8),2);
[Vkr,Hess] = createKrop(@(psi) Hpsi(K,V,psi), v0, 8);
sort(eig(Hess(1:8,:)))
[Pkr,Dkr] = eig(Hess);
[Pkr,Dkr] = eig(Hess(1:8,:));
save