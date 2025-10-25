Htls = [0 0; 0 1]
Hharmonic = diag(0:3)
H0TLSs = kron(Htls, eye(2)) + kron(eye(2), Htls)
kron(Htls, eye(2))
kron(eye(2), Htls)
sigmap = [0 0; 1 0]
coupling = 0.1
sigmam = [0 1; 0 0]
Hcoupling = [kron(sigmap, sigmam) + kron(sigmam, sigmap)]
kron(sigmap, sigmam)
diag([0 1; 1 0])
eig([0 1; 1 0])
[P, D] = eig([0 1; 1 0])
HTLSs = H0TLSs + coupling*Hcoupling
[P, D] = eig(HTLSs)
syms sigmax theta
exp(-1i*theta*[0 1; 1 0])
whos
clear sigmax theta
syms t
exp(-1i*[0 1; 1 0]*t)
exp(-1i*[0 1; 1 0])
sqnorm(ans)
[0 1; 1 0]
[Px, Dx] = eig([0 1; 1 0])
Px*exp(-1i*Dx)*Px'
sqnorm(ans)
exp(-1i*Dx)
Px*expm(-1i*Dx)*Px'
sqnorm(ans)
Px*diag(exp(-1i*diag(Dx)))*Px'
sqnorm(ans)
expm(-1i*[0 1; 1 0]*t)
simplify(ans)
syms omega
expm(-1i*omega/2*[0 1; 1 0]*t)
simplify(ans)
expm(-1i*omega*[0 1; 1 0]*t)
simplify(ans)
whos
UsymTLSs = simplify(expm(-1i*omega*[0 1; 1 0]*t))
save HTLSs
psi0 = [0;1;0;0]
psi = fMdiag(HTLSs, @(H) exp(-1i*H), psi0, 2*pi*20, 1e3);
figure
plot(conj(psi).*psi)
size(psi)
plot(0:40*pi/1e3:40*pi, conj(psi).*psi)
P*H0TLSs*P'
H0TLSs
H0harc = kron(Hharmonic, H0TLSs)
H0harcs = sparse(kron(Hharmonic, H0TLSs))
H0harc = kron(Hharmonic, eye(4)) + kron(eye(4), H0TLSs)
H0harcs = sparse(kron(Hharmonic, eye(4)) + kron(eye(4), H0TLSs))
whos
doc diag
a = diag(sqrt(1:3), 1)
adag = diag(sqrt(1:3), -1)
Htls1 = kron(Htls, eye(2))
Htls2 = kron(eye(2),Htls)
sigmap
sigmap1 = kron(sigmap, eye(2))
sigmam1 = kron(sigmam, eye(2))
sigmap2 = kron(eye(2), sigmap)
sigmam2 = kron(eye(2), sigmam)
Hhcoupling = coupling(kron(adag, sigmam1) + kron(a, sigmap1) + kron(adag, sigmam2) + kron(a, sigmap2))
Hhcoupling = coupling*(kron(adag, sigmam1) + kron(a, sigmap1) + kron(adag, sigmam2) + kron(a, sigmap2))
kron(adag, sigmam1)
kron(a, sigmap1)
Hhcouplings = sparse(coupling*(kron(adag, sigmam1) + kron(a, sigmap1) + kron(adag, sigmam2) + kron(a, sigmap2)))
max(abs(Hhcouplings-Hhcouplings')
max(abs(Hhcouplings-Hhcouplings'))
v=rand(16,1);
Hharcs = H0harcs + Hhcouplings
Hharc = H0harc + Hhcoupling
tic, for j=1:1e3, Hharcs*v; end, toc
tic, for j=1:1e3, Hharc*v; end, toc
tic, for j=1:1e4, Hharc*v; end, toc
tic, for j=1:1e4, Hharcs*v; end, toc
tic, for j=1:1e5, Hharcs*v; end, toc
tic, for j=1:1e5, Hharc*v; end, toc
[Ph, Dh] = eig(Hharcs)
[Ph, Dh] = eig(Hharc)
Eharc = diag(Dh)
Ph*H0harc*Ph'
H0harcfb = Ph*H0harc*Ph';
H0harcfb(abs(H0harcfb)<1e-14)=0
Hharmonic = diag((0:3)*1.1)
Hharmonic = diag(0:3)
Hharmonic1 = diag((0:3)*1.1)
Hharmonic2 = diag((0:3)*2)
H0harc1 = kron(Hharmonic1, eye(4)) + kron(eye(4), H0TLSs);
Hharc1 = H0harc1 + Hhcoupling;
[Ph1, Dh1] = eig(Hharc1);
H0harc1fb = Ph1*H0harc1*Ph1';
H0harc1fb
H0harcfb1(abs(H0harcfb1)<1e-14)=0
H0harc1fb(abs(H0harc1fb)<1e-14)=0
nnz(H0harc1fb)
nnz(H0harcfb)
H0harcfbs = sparse(H0harcfb)
H0harc1fbs = sparse(H0harc1fb)
doc find
H0harcfbs(H0harc1fbs~=0)
H0harcfb(H0harc1fb~=0)
H0harcfbs~=0 + H0harc1fbs~=0
(H0harcfbs~=0) + (H0harc1fbs~=0)
(H0harcfbs==0) & (H0harc1fbs~=0)
H0harc1fbs((H0harcfbs==0) & (H0harc1fbs~=0))
12.5+150-49
H0harcfb
H0harc2 = kron(Hharmonic2, eye(4)) + kron(eye(4), H0TLSs);
Hharc2 = H0harc2 + Hhcoupling;
[Ph2, Dh2] = eig(Hharc2);
H0harc2fb = Ph2*H0harc2*Ph2'
H0harc2fb(abs(H0harc2fb)<1e-14)=0
psih2 = fMdiag(Hharc, @(H) exp(-1i*H), , 2*pi*200, 1e4);
psi0h = zeros(16,1);
psi0h(2) = 1
psih2 = fMdiag(Hharc2, @(H) exp(-1i*H), psi0h2, 2*pi*200, 1e3);
psih2 = fMdiag(Hharc2, @(H) exp(-1i*H), psi0h, 2*pi*200, 1e3);
figure
plot(0:400*pi/1e3:400*pi, conj(psih2).*psih2)
max(conj(psih2).*psih2, 2)
doc max
max(conj(psih2).*psih2,[], 2)
Htls1h = kron(eye(4), Htls1)
Htls2h = kron(eye(4), Htls2)
Htls1hfb = Ph*Htls1h*Ph'
Htls1hfb(Htls1hfb<1e-14) = 0
Htls1hfb = Ph*Htls1h*Ph'
Htls1hfb(abs(Htls1hfb)<1e-14) = 0
Htls2hfb = Ph*Htls2h*Ph'
Htls2hfb(abs(Htls2hfb)<1e-14) = 0
Htlstotfb =  Htls2hfb+Htls1hfb
Htlstotfb(abs(Htlstotfb)<1e-14) = 0
nnz(Htlstotfb)
nnz(Htls1fb)
nnz(Htls1hfb)
nnz(Htls2hfb)
Htls1hfbs = sparse(Htls1hfb)
Htls2hfbs = sparse(Htls2hfb)
(Htls2hfb==0) & (Htls1hfb~=0)
(Htls2hfb==0)
(Htls2hfb~=0)
(Htls2hfb==0) & (Htls1hfb~=0)
Ph
Htls2hfb1 = Ph1*Htls2h*Ph1'
Htls2hfb1(abs(Htls2hfb1)<1e-14) = 0
nnz(Htls2hfb1)
nnz(Htls2hfb)
Htls2hfb1s = sparse(Htls2hfb1)
[Htls2hfbs Htls2hfb1s
]
Htls2hfbs
Htls2hfb1s
abs(Htls2hfbs)>abs(Htls2hfb1s)
Htls1hfb1 = Ph1*Htls1h*Ph1'
Htls1hfb1(abs(Htls1hfb1)>1e-14)=0
Htls1hfb1 = Ph1*Htls1h*Ph1';
Htls1hfb1(abs(Htls1hfb1)<1e-14)=0
nnz(Htls1hfb1)
nnz(Htls2hfb1)
nnz(Htls1hfb1+Htls2hfb1)
(Htls1hfb1+Htls2hfb1)
psih = fMdiag(Hharc, @(H) exp(-1i*H), psi0h, 2*pi*20, 1e3);
figure
plot(0:40*pi/1e3:40*pi, conj(psih).*psih)
max(conj(psih).*psih,[], 2)
psih1 = fMdiag(Hharc1, @(H) exp(-1i*H), psi0h, 2*pi*20, 1e3);
figure
plot(0:40*pi/1e3:40*pi, conj(psih1).*psih1)
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget, H0, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, T, 0.05, 7, 7, 1e-4, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-47.5))).*sin(t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 0:0.1:50, 0.25, 7, 7, 1e-4, 1e3);
target = zeros(16,1);
target(3) = 1
psi0h
t=0:0.1:50;
figure
plot(t, 0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))))
plot(t(1:100), 0.5*(tanh(2*(t(1:100)-2.5))-tanh(2*(t(1:100)-7.5))))
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 1e3*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 0:0.1:50, 0.25, 7, 7, 1e-4, 1e3);
D
Dh
Htls1h
Htls1hs = sparse(Htls1h);
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 1e3*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 0:0.1:50, 0.25, 7, 7, 1e-4, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 1e3*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 50, 0.25, 7, 7, 1e-4, 1e3);
J1
1-J1
figure
plot(0:0.25:50, field)
mean(field)
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 1e2*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 50, 0.25, 7, 7, 1e-4, 1e3);
plot(0:0.25:50, field)
1-J1
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 50*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 50, 0.25, 7, 7, 1e-4, 1e3);
1-J1
figure
plot(0:0.25:50, field)
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 10*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 50, 0.25, 7, 7, 1e-4, 1e3);
1-J1
figure
plot(0:0.25:50, field)
whos
clear invHess psih psih1 psih2
edit evH
figure
plot(0:0.25:50, conj(psi).*psi)
mE = evH(psi, @(psi) Hharcs*psi);
figure
plot(0:0.25:50, mE)
max(conj(psi).*psi,[], 2)
psih1 = fMdiag(Hharc1, @(H) exp(-1i*H), psi0h, 2*pi*20, 1e3);
max(conj(psih1).*psih1,[], 2)
Htls1h
Hharc
eig(Hcoupling)
Hcoupling
eig(Hharc)
whos
Hharc2
superposition = zeros(16,1);
superposition([3 5]) = sqrt(2)/2
superposition'*Hharc1*superposition
superposition'*Hharc*superposition
superposition'*superposition
superposition'*Hharc2*superposition
psih = fMdiag(Hharc, @(H) exp(-1i*H), psi0h, 2*pi*20, 1e3);
max(conj(psih).*psih,[], 2)
figure
plot(0:40*pi/1e3:40*pi, conj(psih).*psih)
diag(Dh)
Ph
Ph(:,2)'*Hharc*Ph(:,2)
Ph(:,3)'*Hharc*Ph(:,3)
Ph(:,4)'*Hharc*Ph(:,4)
2-ans
superposition1 = (Ph(:,2) + Ph(:,4))/sqrt(2)
superposition2 = (Ph(:,2) - Ph(:,4))/sqrt(2)
superposition1'*Hharc*superposition1
superposition2'*Hharc*superposition2
superposition3 = (Ph(:,2) + Ph(:,4))/(2*sqrt(2)) + Ph(:,3)/sqrt(2)
superposition3 = (Ph(:,2) + Ph(:,4))/(2) + Ph(:,3)/sqrt(2)
superposition3'*Hharc*superposition3
Htls1hs
whos
Htls1hfbs
Htls1hfb(1:5, 1:5)
psifb = Ph'*psi;
figure
plot(0:0.25:50, conj(psifb).*psifb)
psifb(:,1)
supersosition3
superposition3
superposition4 = (Ph(:,2) - Ph(:,4))/(2) + Ph(:,3)/sqrt(2)
superposition1
superposition2
psi01 = zeros(16,1);
psi01(2) = 1;
Ph'*psi01
psi10 = zeros(16,1);
psi10(3) = 1;
Ph'*psi10
whos
Eharc
2*pi/(Eharc(3)-Eharc(2))
2*pi/(35.75-14.5)
ans/2
(Eharc(3)-Eharc(2))
mean(field)
Htls1hfb(2:4, 2:4)
max(conj(psifb).*psifb,[], 2)
Ph(:,5)
Htls2h
Htls2hfb
Ph(:,2:4)
Htls1hfb = Ph'*Htls1h*Ph
Htls2hfb = Ph'*Htls1h*Ph
Htls2hfb = Ph'*Htls2h*Ph
Htls1hfb(abs(Htls1hfb)<1e-14) = 0
Htls2hfb(abs(Htls2hfb)<1e-14) = 0
Htls1hfbs = sparse(Htls1hfb)
nnz(Htls1hfbs)
nnz(Htls2hfbs)
Htls1hfb1 = Ph1'*Htls1h*Ph1;
Htls1hfb1(abs(Htls1hfb1)<1e-14)=0
nnz(Htls1hfbs1)
nnz(Htls1hfb1)
Htls2hfbs = sparse(Htls2hfb)
nnz(Htls2hfbs)
whos
Htls2hfb1 = Ph1'*Htls2h*Ph1
Htls2hfb1(abs(Htls2hfb1)<1e-14) = 0
Htls2hfb1s = sparse(Htls2hfb1)
Htlstotfb =  Htls2hfb+Htls1hfb
Htlstotfb(abs(Htlstotfb)<1e-14) = 0
save HTLSs_harmonic
Hharmonic3 = diag((0:3)*0.9)
H0harc2fb = Ph2'*H0harc2*Ph2
H0harc2fb(abs(H0harc2fb)<1e-14)=0
H0harc3 = kron(Hharmonic3, eye(4)) + kron(eye(4), H0TLSs);
Hharc3 = H0harc3 + Hhcoupling;
Hhcoupling
[Ph3, Dh3] = eig(Hharc3);
Eharc3 = diag(Dh3)
2*pi/(35.75-14.5)
%-- 05/08/2020 10:36 --%
load HTLSs_harmonic
%-- 05/08/2020 13:03 --%
load HTLSs_harmonic
Eharc1 = diag(Dh1)
Htls1hfb1(2:4, 2:4)
Hharmonic3 = diag((0:3)*0.9)
H0harc2fb = Ph2'*H0harc2*Ph2
H0harc2fb(abs(H0harc2fb)<1e-14)=0
H0harc3 = kron(Hharmonic3, eye(4)) + kron(eye(4), H0TLSs);
Hharc3 = H0harc3 + Hhcoupling;
Hhcoupling
[Ph3, Dh3] = eig(Hharc3);
Eharc3 = diag(Dh3)
2*pi/(35.75-14.5)
save HTLSs_harmonic
fieldw = dctI(field);
figure
plot(0:pi/50:pi, fieldw(1:51))
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))).*sin(t), @(t) 1e3*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-37.5))), options, 50, 0.25, 7, 7, 1e-4, 1e3);
1-J11
figure
plot(0:0.25:50, field1)
2*pi/3
fieldw1 = dctI(field1);
figure
plot(0:pi/50:pi, fieldw1(1:51))
[allfield2, field2, psi2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*sin(t), @(t) 20*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
J12
1-J12
figure
plot(0:0.25:100, field2)
figure
fieldw2 = dctI(field2);
plot(0:pi/100:pi, fieldw2(1:101))
Hharc
coupling
psi = fMdiag(HTLSs, @(H) exp(-1i*H), psi0, 2*pi*20, 1e3);
figure
plot(0:40*pi/1e3:40*pi, conj(psi).*psi)
figure
plot(0:0.25:50, field)
1-J1
1-J12
UsymTLSs
[allfield3, field3, psi3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, alpha3, invHess3] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 20*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
1-J13
figure
plot(0:0.25:70, field3)
[allfield3, field3, psi3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, alpha3, invHess3] = OCqnTDpenal(psi0h, target, Hharcs, [-1 6], Htls1hs, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
1-J13
plot(0:0.25:70, field3)
fieldw3 = dctI(field3);
figure
plot(0:pi/70:pi, fieldw3(1:101))
plot(0:pi/70:pi, fieldw3(1:71))
2*pi/20
Hharc
Hharc([2,3,5],[2,3,5])
Htls1h([2,3,5],[2,3,5])
[allfield3a, field3a, psi3a, relE3a, conv3a, niter3a, mallniterc3a, J13a, maxgrad3a, alpha3a, invHess3a] = OCqnTDpenal(psi0h([2,3,5]), target([2,3,5]), Hharcs([2,3,5],[2,3,5]), [-1 6], Htls1hs([2,3,5],[2,3,5]), @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
max(abs(field3-field3a))
niter3
psi3([2,3,5],end) - psi3a(:,end)
psi3a
psi3a(:,end)
psi3a(:,end).*conj(psi(:,end))
psi3a(:,end).*conj(psi3a(:,end))
Htls1h3 = Htls1h([2,3,5],[2,3,5]);
Hharc3 = Hharc([2,3,5],[2,3,5])
target3 = target([2,3,5])
psi0h3  = psi0h([2,3,5])
psih3u =zeros(6,1);
psih3u([1,5]) = 1
target3u =zeros(6,1);
target3u([2,4]) = 1
Hharc3u = kron(eye(2), Hharc3)
Htls1h3u = kron(eye(2), Htls1h3)
[allfieldu3, fieldu3, psiu3, relEu3, convu3, niteru3, mallnitercu3, J1u3, maxgradu3, alphau3, invHessu3] = OCqnTDpenal(psi0h3u, target3u, Hharcs3u, [-1 3], Htls1hs3u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
psih03u = psi3u
psih03u = psih3u
clear(psih3u)
clear psih3u
[allfieldu3, fieldu3, psiu3, relEu3, convu3, niteru3, mallnitercu3, J1u3, maxgradu3, alphau3, invHessu3] = OCqnTDpenal(psi0h3u, target3u, Hharcs3u, [-1 3], Htls1hs3u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu3, fieldu3, psiu3, relEu3, convu3, niteru3, mallnitercu3, J1u3, maxgradu3, alphau3, invHessu3] = OCqnTDpenal(psih03u, target3u, Hharcs3u, [-1 3], Htls1hs3u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu3, fieldu3, psiu3, relEu3, convu3, niteru3, mallnitercu3, J1u3, maxgradu3, alphau3, invHessu3] = OCqnTDpenal(psih03u, target3u, Hharc3u, [-1 3], Htls1h3u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
4-J1u3
1-J13
2*ans
figure
plot(0:0.25:70, fieldu3)
max(abs(fieldu3-fieldu3(end:-1:1)))
figure
plot(0:0.25:70, psiu3.*conj(psiu3))
plot(0:0.25:70, psiu3(1:3,:).*conj(psiu3(1:3,:)))
figure
plot(0:0.25:70, psiu3(4:6,:).*conj(psiu3(4:6,:)))
2*pi/17
Hharc
Hharcfb
H0harcfb
whos
fieldwu3 = dctI(fieldu3);
figure
plot(0:pi/70:pi, fieldwu3(1:71))
hold on
plot(0:pi/70:pi, fieldw3(1:71))
Htls1hfb
Ph1(:, 1:9)
Hharc8 = Hharc([1:7,9],[1:7,9])
Htls1h3 = Htls1h([1:7,9],[1:7,9]);
Htls1h3 = Htls1h([1:7,9],[1:7,9])
Htls1h3 = Htls1h([2,3,5],[2,3,5]);
Htls1h8 = Htls1h([1:7,9],[1:7,9])
psih08u = [eye(4);zeros(4)];
psih08u(2:3,2:3) = [0 1i; 1i 0]
psih08u = [eye(4);zeros(4)];
target8u = psih08u;
target8u(2:3,2:3) = [0 1i; 1i 0]
Hharc8u = kron(sparse(eye(4)), sparse(Hharc8))
Htls1h8u = kron(sparse(eye(4)), Htls1h3)
Htls1h8u = kron(sparse(eye(4)), sparse(Htls1h8))
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-1 3], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8
figure
plot(0:0.25:70, fieldu8)
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8
figure
plot(0:0.25:70, fieldu8)
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(0.14t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
Eharc
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*sin(0.1414*t), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*cos(0.1414*(t-35)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8
figure
plot(0:0.25:70, fieldu8)
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, allfieldu3, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 70, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8
figure
plot(0:0.25:70, fieldu8)
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, allfieldu3, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))).*cos(0.1414*(t-35)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-57.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8
figure
plot(0:0.25:100, fieldu8)
psiu8(:,end).*conj(psiu8(:,end))
reshape(psiu8(:,end).*conj(psiu8(:,end)), 8, 4)
reshape(psiu8(:,end), 8, 4)
target8u = psih08u;
target8u(2:3,2:3) = [0 1i; 1i 0]
target8u1 = psih08u;
target8u1(2:3,2:3) = [0 1; 1 0]
[allfieldu8a, fieldu8a, psiu8a, relEu8a, convu8a, niteru8a, mallnitercu8a, J1u8a, maxgradu8a, alphau8a, invHessu8a] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8a
figure
plot(0:0.25:100, fieldu8a)
options1 = optionsOCqn(1e-6, 1e3);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
options1.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
options.invHess0 = invHessu8a;
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfield8ua, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-4, 1e3);
options.invHess0 = [];
options1.invHess0 = invHessu8a;
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-4, 1e3);
options1.invHess0 = [];
options1.invHess0 = invHessu8a;
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-6, 1e3);
options1 = optionsOCqn(1e-8, 1e3);
options1.invHess0 = invHessu8a;
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-8, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.1, 7, 7, 1e-8, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 9, 9, 1e-8, 1e3);
J1u8a1
16-ans
16-J1u8a
16-J1u8a1
hold on
plot(0:0.25:100, fieldu8a1)
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.01*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
figutr
figure
plot(0:0.25:100, fieldu8a1)
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 99*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 9, 9, 1e-8, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 99*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 90*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 99*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-8, 1e3);
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 90*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options1, 100, 0.25, 7, 7, 1e-8, 1e3);
16-J1u8a
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 10*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8a1
figure
plot(0:0.25:100, fieldu8a1)
hold on
plot(0:0.25:100, fieldu8a)
[allfieldu8a1, fieldu8a1, psiu8a1, relEu8a1, convu8a1, niteru8a1, mallnitercu8a1, J1u8a1, maxgradu8a1, alphau8a1, invHessu8a1] = OCqnTDpenal(psih08u(:), target8u1(:), Hharc8u, [-1 3], Htls1h8u, allfieldu8a, @(t) 10*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8a1
16-J1u8
figure
plot(0:0.25:100, fieldu8)
[allfieldu8b, fieldu8b, psiu8b, relEu8b, convu8b, niteru8b, mallnitercu8b, J1u8b, maxgradu8b, alphau8b, invHessu8b] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 10*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8b
[allfieldu8b, fieldu8b, psiu8b, relEu8b, convu8b, niteru8b, mallnitercu8b, J1u8b, maxgradu8b, alphau8b, invHessu8b] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))).*cos(0.1414*(t-100)), @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options, 200, 0.25, 7, 7, 1e-4, 1e3);
16-J1u8b
figure
plot(0:0.25:100, fieldu8b)
plot(0:0.25:200, fieldu8b)
options2 = optionsOCqn(1e-6, 2e3);
options2.invHess0 = invHessu8b;
options2 = optionsOCqn(1e-6, 1e4);
options2.invHess0 = invHessu8b;
[allfieldu8b1, fieldu8b1, psiu8b1, relEu8b1, convu8b1, niteru8b1, mallnitercu8b1, J1u8b1, maxgradu8b1, alphau8b1, invHessu8b1] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, allfieldu8b, @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options2, 200, 0.25, 7, 7, 1e-4, 1e3);
[allfieldu8b1, fieldu8b1, psiu8b1, relEu8b1, convu8b1, niteru8b1, mallnitercu8b1, J1u8b1, maxgradu8b1, alphau8b1, invHessu8b1] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, allfieldu8b, @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options2, 200, 0.25, 7, 7, 1e-6, 1e4);
figure
plot(0:0.25:200, fieldu8b1)
16-J1u8b
16-J1u8b1
figure
200/(5*2*pi)
plot(0:0.25:200, conj(psiu8b1).*psiu8b1)
plot(0:0.25:200, conj(psiu8b1(1:8)).*psiu8b1(1:8))
plot(0:0.25:200, conj(psiu8b1(1:8,:)).*psiu8b1(1:8,:))
plot(0:0.25:200, conj(psiu8b1(9:16,:)).*psiu8b1(9:16,:))
figure
plot(0:0.25:200, conj(psiu8b1(17:24,:)).*psiu8b1(17:24,:))
figure
plot(0:0.25:200, conj(psiu8b1(25:32,:)).*psiu8b1(25:32,:))
whos
clear invHess1 invHess2 invHess3 invHess3a invHessu3 invHessu8 invHessu8a invHessu8a1 invHessu8b invHessu8b1 options2
whos
clear options1
whos
save HTLSs_harmonic
H0
whos
HTLSs
H03 = diag([0 1 2-(0.2/5)])
0.2/5
Hcoupling
a
adag
%-- 24/08/2020 10:52 --%
load HTLSs_harmonic
H03 = diag([0 1 2-(0.2/5)])
H03LSs = kron(H03, eye(3)) + kron(eye(3), H03)
Hharmonic
H0harc = kron(Hharmonic, eye(9)) + kron(eye(9), H03LSs)
H0harc = kron(Hharmonic, eye(9)) + kron(eye(4), H03LSs)
H0harcs = sparse(kron(Hharmonic, eye(9)) + kron(eye(4), H03LSs))
a
a3 = a(1:3,1:3)
adag3 = adag(1:3,1:3)
sigmap31 = kron(adag3, eye(3))
sigmam31 = kron(a3, eye(3))
sigmap32 = kron(eye(3), adag3)
sigmam32 = kron(eye(3), a3)
Hhcoupling = coupling*(kron(adag, sigmam1) + kron(a, sigmap1) + kron(adag, sigmam2) + kron(a, sigmap2))
H3hcoupling = coupling*(kron(adag, sigmam31) + kron(a, sigmap31) + kron(adag, sigmam32) + kron(a, sigmap32))
H3hcouplings = sparse(coupling*(kron(adag, sigmam31) + kron(a, sigmap31) + kron(adag, sigmam32) + kron(a, sigmap32)))
H3hcouplings-H3hcouplings'
whos
H0harc = kron(Hharmonic, eye(4)) + kron(eye(4), H0TLSs);
H0harcs = sparse(kron(Hharmonic, eye(4)) + kron(eye(4), H0TLSs));
H03harc = kron(Hharmonic, eye(9)) + kron(eye(4), H03LSs);
H03harcs = sparse(kron(Hharmonic, eye(9)) + kron(eye(4), H03LSs));
save HTLSs_harmonic
%-- 31/08/2020 12:27 --%
load HTLSs_harmonic
Hharc = H03harc + H3hcoupling
Hharcs = H03harcs + H3hcouplings
[Ph, Dh] = eig(Hharc)
[Ph3, Dh3] = eig(H3harc);
Hharc = H0harc + Hhcoupling
Hharcs = H0harcs + Hhcouplings;
H3harcs = H03harcs + H3hcouplings;
H3harc = H03harc + H3hcoupling;
H3harcs
[Ph3, Dh3] = eig(H3harc);
E3harc = diag(Dh3);
E3harc = diag(Dh3)
psi3h = fMdiag(Hharc3, @(H) exp(-1i*H), , 2*pi*20, 1e3);
psi0h = zeros(36,1);
psi0h = zeros(16,1);
psi03h = zeros(36,1);
psi03h(2) = 1;
psi3h = fMdiag(Hharc3, @(H) exp(-1i*H), psi03h, 2*pi*20, 1e3);
psi3h = fMdiag(H3harc, @(H) exp(-1i*H), psi03h, 2*pi*20, 1e3);
figure
plot(0:40*pi/1e3:40*pi, conj(psi3h).*psi3h)
conj(psi3h(:,end)).*psi3h(:, end)
psi03ha = zeros(36,1);
psi03ha(5) = 1;
psi3ha = fMdiag(H3harc, @(H) exp(-1i*H), psi03ha, 2*pi*20, 1e3);
figure
plot(0:40*pi/1e3:40*pi, conj(psi3ha).*psi3ha)
conj(psi3ha(:,end)).*psi3ha(:, end)
psi3ha_f = ans;
psi3ha_f(psi3ha_f<1e-14) = 0;
sparse(psi3ha_f)
H3harc
isingle3 = [2 4 10];
idouble3 = [3 5 7 11 13 19];
save HTLSs_harmonic
%-- 01/09/2020 7:23 --%
load HTLSs_harmonic
H3LS_1 = kron(H03, eye(3))
H3LS_2 = kron(eye(3), H03)
H3LS_1h = kron(eye(4), H3LS_1);
H3LS_1hs = sparse(kron(eye(4), H3LS_1));
H3LS_1hs = sparse(kron(eye(4), H3LS_1))
H3harcu = zeros(13);
H3harcu(2:4, 2:4) = H3harc(isingle3, isingle3);
H3harcu(5:7, 5:7) = H3harc(isingle3, isingle3);
H3harcu(8:13, 8:13) = H3harc(idouble3, idouble3);
H3harcsu
H3harcu
H3harcus = sparse(H3harcu);
H3LS_1hu = zeros(13);
clear options1
whos
save HTLSs_harmonic
H3LS_1hu = zeros(13);
H3LS_1hu(2:4, 2:4) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(5:7, 5:7) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(8:13, 8:13) = H3LS_1h(idouble3, idouble3)
target3LSu = zeros(13,1);
target3LSu(1) = 1;
target3LSu(3) = 1i;
target3LSu(5) = 1i;
target3LSu(9) = 1
u0_3 = zeros(13);
u0_3([1 2 6 9])=1
u0_3 = zeros(13,1);
u0_3([1 2 6 9])=1
[u0_3, target3LSu]
[allfieldu8, fieldu8, psiu8, relEu8, convu8, niteru8, mallnitercu8, J1u8, maxgradu8, alphau8, invHessu8] = OCqnTDpenal(psih08u(:), target8u(:), Hharc8u, [-2 6], Htls1h8u, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
H3LS_1hus = sparse(H3LS_1hus)
H3LS_1hus = sparse(H3LS_1hu)
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.1, 7, 7, 1e-4, 1e3);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 9, 9, 1e-4, 1e3);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
options2 = options;
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 5);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
options2.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 5);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.25, 7, 7, 1e-4, 1e3);
options2.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 3);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options2, 100, 0.25, 7, 7, 1e-4, 1e3);
options2.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 5);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options2, 100, 0.25, 7, 7, 1e-4, 1e3);
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 100*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))), options, 100, 0.1, 7, 7, 1e-4, 1e3);
u3(:, end)
u3(:, end).*conj(u3(:,end))
16-J13u
figure
plot(0:0.1:100, field3u)
[allfield3u, field3u, u3, relE3u, conv3u, niter3u, mallniterc3u, J13u, maxgrad3u, alpha3u, invHess3u] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-87.5))).*cos(0.1414*(t-50)), @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options, 200, 0.25, 7, 7, 1e-4, 1e3);
16-J13u
options2 = options;
options2.invHess0 = invHess3u;
[allfield3ua, field3ua, u3a, relE3ua, conv3ua, niter3ua, mallniterc3ua, J13ua, maxgrad3ua, alpha3ua, invHess3ua] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, allfield3u, @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options2, 200, 0.25, 7, 7, 1e-4, 1e3);
options3 = options;
options3.invHess0 = invHess3ua;
options3 = optionsOCqn(1e-4, 1e4);
16-J13ua
plot(0:0.25:200, field3u)
plot(0:0.25:200, field3ua)
[allfield3ub, field3ua, u3a, relE3ua, conv3ua, niter3ua, mallniterc3ua, J13ua, maxgrad3ua, alpha3ua, invHess3ua] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, allfield3ua, @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options2, 200, 0.25, 7, 7, 1e-4, 1e3);
options3.invHess0 = invHess3ua;
[allfield3ub, field3ub, u3b, relE3ub, conv3ub, niter3ub, mallniterc3ub, J13ub, maxgrad3ub, alpha3ub, invHess3ub] = OCqnTDpenal(u0_3, target3LSu(:), H3harcus, [-2 6], H3LS_1hus, allfield3ua, @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options3, 200, 0.25, 7, 7, 1e-4, 1e3);
16-J13ub
16-J13ua
figure
plot(0:0.25:200, field3ub)
figure
plot(0:0.25:200, u3b.*conj(u3b))
plot(0:0.25:200, u3b(1,:).*conj(u3b(1,:)))
plot(0:0.25:200, u3b(2:4,:).*conj(u3b(2:4,:)))
plot(0:0.25:200, u3b(5:7,:).*conj(u3b(5:7,:)))
plot(0:0.25:200, u3b(8:13,:).*conj(u3b(8:13,:)))
u3b(:,end).*conj(u3b(:,end))
u3b(:,end)
max(conj(u3b).*u3b, 2)
max(conj(u3b).*u3b, [], 2)
whos
clear invHess3u invHess3ua invHess3ub options2 options3
whos
save HTLSs_harmonic
options4 = optionsOCqn(1e-6, 1e4);
clear options4
size(allfield)
coupling
0.15/5
coupling1 = 0.15/5
0.4/5
0.3/5
options4 = optionsOCqn(1e-6, 1e4);
H3hcoupling1 = coupling1*(kron(adag, sigmam31) + kron(a, sigmap31) + kron(adag, sigmam32) + kron(a, sigmap32))
H3harc1 = H03harc + H3hcoupling1;
H3harc1s = sparse(H3harc1)
[allfield3u1, field3u1, u31, relE3u1, conv3u1, niter3u1, mallniterc3u1, J13u1, maxgrad3u1, alpha3u1, invHess3u1] = OCqnTDpenal(u0_3, target3LSu(:), H3harc1us, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))).*cos(0.1414*(t-50)), @(t) 200*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-187.5))), options, 200, 0.25, 7, 7, 1e-4, 1e3);
H3harc1u = zeros(13);
H3harc1u(2:4, 2:4) = H3harc1(isingle3, isingle3);
H3harc1u(5:7, 5:7) = H3harc1(isingle3, isingle3);
H3harc1u(8:13, 8:13) = H3harc1(idouble3, idouble3)
H3harc1us = sparse(H3harc1u)
3*sqrt(2)
[allfield3u1, field3u1, u31, relE3u1, conv3u1, niter3u1, mallniterc3u1, J13u1, maxgrad3u1, alpha3u1, invHess3u1] = OCqnTDpenal(u0_3, target3LSu(:), H3harc1us, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-387.5))).*cos(0.1414*(t-50)), @(t) 400*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-387.5))), options4, 400, 0.25, 7, 7, 1e-4, 1e3);
[allfield3u1, field3u1, u31, relE3u1, conv3u1, niter3u1, mallniterc3u1, J13u1, maxgrad3u1, alpha3u1, invHess3u1] = OCqnTDpenal(u0_3, target3LSu(:), H3harc1us, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-387.5))).*cos(0.1414*(t-50)), @(t) 400*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-787.5))), options4, 800, 0.25, 7, 7, 1e-4, 1e3);
[allfield3u1, field3u1, u31, relE3u1, conv3u1, niter3u1, mallniterc3u1, J13u1, maxgrad3u1, alpha3u1, invHess3u1] = OCqnTDpenal(u0_3, target3LSu(:), H3harc1us, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-787.5))).*cos(0.1414*(t-400)), @(t) 400*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-787.5))), options4, 800, 0.25, 7, 7, 1e-4, 1e3);
options4.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
[allfield3u1, field3u1, u31, relE3u1, conv3u1, niter3u1, mallniterc3u1, J13u1, maxgrad3u1, alpha3u1, invHess3u1] = OCqnTDpenal(u0_3, target3LSu(:), H3harc1us, [-2 6], H3LS_1hus, @(t) 0.1*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-787.5))).*cos(0.1414*(t-400)), @(t) 400*0.5*(tanh(0.4*(t-12.5))-tanh(0.4*(t-787.5))), options4, 800, 0.25, 7, 7, 1e-4, 1e3);
relE3u1
16-J13u1
max(abs(allfield3u1))
figure
plot(0:0.25:800, field3u1)
Hcoupling1
coupling1
%-- 15/09/2020 17:37 --%
load HTLSs_harmonic
M3 = zeros(2,2,3)
v=ones(3,1);
v*M3
v.*M3
v(1:2)*M3(:,:,1)
v(1:2).*M3(:,:,1)
M=[1 2; 3 4]
v(1:2).*M3(:,:,1)
V=[5;6]
v.*M
v=[5;6]
v.*M
v.'.*M
v.*v
v.*v.'
clear all
psi0 = [1;0;0];
target = [0;1;0];
H0 = diag(2*pi*[0, 5, 9.8]);
Edomain = 2*pi*[-5, 15];
miu = [0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0];
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
[Hoperations, fcouplingOp] = Hmats2Hops(H0, miu)
ftarget = phi2ftarget(target)
fguess
whos
figure
plot(0:0.01:10, fieldu1)
fieldu1w = dctI(fieldu1w);
fieldu1w = dctI(fieldu1)*10/(sqrt(1e3*pi));
figure
plot(0:pi/10:pi/0.01, fieldu1w)
plot(0:pi/10:20*pi, fieldu1w(1:201))
hold on
w=0:pi/10:pi/0.01;
plot(w(1:201), 0.5*(tanh(w(1:201)-20*pi)+1))
plot(w(1:251), 0.5*(tanh(w(1:251)-20*pi)+1))
plot(w(1:251), 0.5*(1-tanh(w(1:251)-20*pi)))
plot(w(1:201), exp((w(1:201)-10*pi).^2/(2*25)))
plot(w(1:201), exp(-(w(1:201)-10*pi).^2/(2*25)))
plot(w(1:201), exp(-(w(1:201)-10*pi).^2/(2*9)))
plot(w(1:201), exp(-(w(1:201)-10*pi).^2/(2*9)).*sin((w(1:201)-10*pi)))
filterE = @(w) 0.5*(1-tanh(w(1:251)-20*pi))
filterE = @(w) 0.5*(1-tanh(w-20*pi))
fguess = @(w) exp(-(w-10*pi).^2/(2*9)).*sin((w-10*pi))
plot(w(1:201), fguess(w(1:201)))
plot(w(1:201), filterE(w(1:201)))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE, options, 10, 0.01, 7, 7, 1e-7);
filterE = @(w) 1e3*0.5*(1-tanh(w-20*pi))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE, options, 10, 0.01, 7, 7, 1e-7);
Jterms
1-Jterms.Jmax
figure
plot(w(1:201), fieldw(1:201))
plot(w, fieldw)
figure
t=0:0.01:10;
plot(t, fieldt)
filterE1 = @(w) 1e3*exp(-(w-10*pi).^2/(2*9)).*sin((w-10*pi))
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
filterE1 = @(w) 1e2*exp(-(w-10*pi).^2/(2*9)).*sin((w-10*pi))
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 20);
filterE1 = @(w) 1e3*exp(-(w-10*pi).^2/(2*9)).*sin((w-10*pi))
options.f_max_alpha = get_f_max_alphaOCf(10, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
edit get_f_max_alphaOCf
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(20, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
filterE1
filterE1 = @(w) 1e2*exp(-(w-10*pi).^2/(2*9)).*sin((w-10*pi))
options.f_max_alpha = get_f_max_alphaOCf(20, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
figure
plot(w(1:201), filterE1(w(1:201)))
filterE1 = @(w) 1e3*exp(-(w-10*pi).^2/(2*9)).*sin((w-10*pi))
fguess
options.f_max_alpha = get_f_max_alphaOCf(100, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
options.f_max_alpha = get_f_max_alphaOCf(10, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
direction
fgrad0
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(10, 0.01, 10, filterE1);
target
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
max(abs(chi_coup_psi))
max(max(abs(allchi)))
(max(abs(chiT)))
max(abs(psiT))
psiT
max(abs(allfield))
max(abs(fieldw))
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
fguess(w(iEnz))
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
fguess(w(wnzi))
fieldw(wnzi)
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
max(abs(fieldw))
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
max(abs(fieldw))
max(abs(fieldwnz))
filterE1 = @(w) 1e3*exp(-(w-10*pi).^2/(2*9))
options.f_max_alpha = get_f_max_alphaOCf(10, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
options.f_max_alpha = get_f_max_alphaOCf(20, 0.01, 10, filterE1);
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
1-Jterms1.Jmax
1-Jterms.Jmax
figure
plot(w(1:201), fieldw1(1:201))
figure
plot(t, fieldt1)
A = [1 2; 3 4];
B= [4 5; 6 7];
fun2 = [@(v) A*v; @(v) B*v]
doc cell
fun2 = {@(v) A*v; @(v) B*v}
whos
size(fun2)
fun2[1]
fun2{1}
fun2{2}
fun2{1}(v)
fun2{1}([3 3])
fun2{1}([3; 3])
v
v=[1;2;3]
repmat(v,3,1)
whos
isa(fun2, function_handle)
isa(fun2, 'function_handle')
isa(fun2{1}, 'function_handle')
isa(fun2, 'cell')
v
v(1:3) = 1
v=[1;2;3]
M=zeros(3);
M(:, :) = v
M(:, :) = repmat(v, 1, 3)
M(:, :) = v.*ones(3)
u = 3:5
v
v(1:3) = u
M
M=M(:)
reshape(M, 3,3)
reshape(M, 2,3)
reshape([M v], 2,8)
reshape([M v], 2,6)
[M v]
M
M= reshape([M v], 3,3)
M= reshape([M], 3,3)
reshape([M v], 2,8)
reshape([M v], 2,6)
ans+[2;3]
[fieldt1a, fieldw1a, psi1a, relE1a, conv1a, niter1a, mallniterc1a, Jterms1a, maxgrad1a, alpha1a, invHess1a] = OClimf_qn(psi0, ftarget, Hoperations, 1, fcouplingOp, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
max(abs(fieldt1-fieldt1a))
lenght(conv)
length(conv)
length(conv1)
whos
miu
P = [0 -1*1i 0; 1i 0 -sqrt(2)*1i; 0 sqrt(2)*1i 0]
P-P'
[Hoperations2, fcouplingOp2] = Hmats2Hops2(H0, miu, P)
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
if [true false] v(3) = 5; end
v
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options, 10, 0.01, 7, 7, 1e-7);
edit dctI
options2 = optionsOCqn(1e-4, 1e3);
options2.f_max_alpha = get_f_max_alphaOCf_multiE(20, 0.01, 10, filterE1);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options2, 10, 0.01, 7, 7, 1e-7);
options2.f_max_alpha = get_f_max_alphaOCf_multiE(20, 0.01, 10, filterE1, 2);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options2, 10, 0.01, 7, 7, 1e-7);
options2.f_max_alpha = get_f_max_alphaOCf_multiE(20, 0.01, 10, filterE1, 2);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options2, 10, 0.01, 7, 7, 1e-7);
options2.f_max_alpha = get_f_max_alphaOCf_multiE(20, 0.01, 10, filterE1, 2);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options2, 10, 0.01, 7, 7, 1e-7);
options2.f_max_alpha = get_f_max_alphaOCf_multiE(20, 0.01, 10, filterE1, 2);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_qn(psi0, ftarget, Hoperations2, 2, fcouplingOp2, Edomain, fguess, filterE1, options2, 10, 0.01, 7, 7, 1e-7);
Jterms2
1-Jterms2
1-Jterms2.Jmax
figure
plot(t, fieldt2)
size(fieldt2)
size(fieldw2)
plot(w(1:201), fieldw2(:, 1:201))
figure
plot(t, dctI(fieldw2(1,:))
plot(t, dctI(fieldw2(1,:)))
hold on
plot(t, dctI(fieldw2(2,:)))
load HTLSs_harmonic
whos
figure
coupling1 = 0.15/5
H3hcoupling1 = coupling1*(kron(adag, sigmam31) + kron(a, sigmap31) + kron(adag, sigmam32) + kron(a, sigmap32))
H3harc1 = H03harc + H3hcoupling1;
H3harc1s = sparse(H3harc1)
H3harc1u = zeros(13);
H3harc1u(2:4, 2:4) = H3harc1(isingle3, isingle3);
H3harc1u(5:7, 5:7) = H3harc1(isingle3, isingle3);
H3harc1u(8:13, 8:13) = H3harc1(idouble3, idouble3)
H3harc1us = sparse(H3harc1u)
options4 = optionsOCqn(1e-6, 1e4);
options4.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
size(field3ub)
load HTLSs_harmonic
size(field3ub)
whos
size(field3u)
size(fieldu3)
H3harc1u = zeros(13);
H3harc1u(2:4, 2:4) = H3harc1(isingle3, isingle3);
H3harc1u(5:7, 5:7) = H3harc1(isingle3, isingle3);
H3harc1u(8:13, 8:13) = H3harc1(idouble3, idouble3)
H3harc1us = sparse(H3harc1u)
eig(H3harc1u)
whos
figure
plot(0:pi/70:pi, fieldw3(1:71))
plot(0:pi/10:20*pi, fieldu1w(1:201))
plot(0:pi/100:pi, fieldw2(1:101))
pi/0.25
options4
options4 = optionsOCqn(1e-4, 1e4);
options4.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterEgates);
filterE = @(w) 1e2*0.5*(1-tanh(w-1))
filterE = @(w) 1e3*0.5*(1-tanh(w-20*pi))
filterE3u1 = @(w)1e2*0.5*(1-tanh(w-1))
w3u1 = 0:pi/800:pi/0.25;
figure
plot(w3u1(1:800), filterE3u1(1:800))
filterE3u1 = @(w)1e2*0.5*(1-tanh(w))
plot(w3u1(1:800), filterE3u1(1:800))
filterE3u1 = @(w)1e2*0.5*(1-tanh(w-1))
plot(w3u1(1:800), filterE3u1(w3u1(1:800)))
filterE3u1 = @(w)1e2*0.5*(1-tanh(5*(w-1)))
plot(w3u1(1:800), filterE3u1(w3u1(1:800)))
filterE3u1 = @(w)1e2*0.5*(1-tanh(10*(w-1)))
plot(w3u1(1:800), filterE3u1(w3u1(1:800)))
options4.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterE3u1);
[Hoperations3u1, fcouplingOp3u1] = Hmats2Hops(H3harc1us, H3LS_1hus);
H3LS_1hu = zeros(13);
H3LS_1hu(2:4, 2:4) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(5:7, 5:7) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(8:13, 8:13) = H3LS_1h(idouble3, idouble3)
target3LSu = zeros(13,1);
H3LS_1 = kron(H03, eye(3))
H3LS_2 = kron(eye(3), H03)
H3LS_1h = kron(eye(4), H3LS_1);
H3LS_1hs = sparse(kron(eye(4), H3LS_1));
H3LS_1hs = sparse(kron(eye(4), H3LS_1))
H3harcu = zeros(13);
H3harcu(2:4, 2:4) = H3harc(isingle3, isingle3);
H3harcu(5:7, 5:7) = H3harc(isingle3, isingle3);
H3harcu(8:13, 8:13) = H3harc(idouble3, idouble3);
H3harcsu
H3harcu
H3harcus = sparse(H3harcu);
H3LS_1hu = zeros(13);
H3harcus = sparse(H3harcu);
H3LS_1hu = zeros(13);
H3LS_1hu(2:4, 2:4) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(5:7, 5:7) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(8:13, 8:13) = H3LS_1h(idouble3, idouble3)
target3LSu = zeros(13,1);
target3LSu(1) = 1;
target3LSu(3) = 1i;
target3LSu(5) = 1i;
target3LSu(9) = 1
u0_3 = zeros(13);
u0_3([1 2 6 9])=1
u0_3 = zeros(13,1);
u0_3([1 2 6 9])=1
[u0_3, target3LSu]
[Hoperations3u1, fcouplingOp3u1] = Hmats2Hops(H3harc1us, H3LS_1hus);
H3LS_1hus = sparse(H3LS_1hus)
H3LS_1hus = sparse(H3LS_1hu)
H3LS_1hus = sparse(H3LS_1hu)
[Hoperations3u1, fcouplingOp3u1] = Hmats2Hops(H3harc1us, H3LS_1hus);
ftarget3u1 = phi2ftarget(target3LSu(:))
[fieldt3u1, fieldw3u1, psi3u1, relE3u1, conv3u1, niter3u1, mallniterc3u1, Jterms3u1, maxgrad3u1, alpha3u1, invHess3u1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1, options4, 800, 0.25, 7, 7, 1e-6)
max(abs(allfield))
Jterms3u1
16-Jterms3u1.Jmax
[fieldt3u1, fieldw3u1, psi3u1, relE3u1, conv3u1, niter3u1, mallniterc3u1, Jterms3u1, maxgrad3u1, alpha3u1, invHess3u1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1, options4, 800, 0.25, 7, 7, 1e-6)
filterE3u1 = @(w)10*0.5*(1-tanh(10*(w-1)))
[fieldt3u1, fieldw3u1, psi3u1, relE3u1, conv3u1, niter3u1, mallniterc3u1, Jterms3u1, maxgrad3u1, alpha3u1, invHess3u1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1, options4, 800, 0.25, 7, 7, 1e-6);
16-Jmax
max(abs(allfield))
16-Jterms3u1.Jmax
figure
plot(0:0.25:800, fieldt)
plot(0:0.25:800, fieldt3u1)
figure
plot(0:pi/800:pi, fieldw3u1(1:801))
filterE3u1a = @(w)10*0.5*(1-tanh(10*(w-0.5)))
figure
plot
plot(w3u1(1:801), filterE3u1a(w3u1(1:801)))
filterE3u1a = @(w)10*0.5*(1-tanh(10*(w-0.25)))
plot(w3u1(1:801), filterE3u1a(w3u1(1:801)))
filterE3u1a = @(w)10*0.5*(1-tanh(20*(w-0.25)))
plot(w3u1(1:801), filterE3u1a(w3u1(1:801)))
figure
plot(0:pi/800:pi, fieldw3u1(1:801))
[fieldt3u1a, fieldw3u1a, psi3u1a, relE3u1a, conv3u1a, niter3u1a, mallniterc3u1a, Jterms3u1a, maxgrad3u1a, alpha3u1a, invHess3u1a] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1a, options4, 800, 0.25, 7, 7, 1e-6);
options4a = optionsOCqn(1e-4, 1e4);
options4a.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterE3u1a);
[fieldt3u1a, fieldw3u1a, psi3u1a, relE3u1a, conv3u1a, niter3u1a, mallniterc3u1a, Jterms3u1a, maxgrad3u1a, alpha3u1a, invHess3u1a] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1a, options4, 800, 0.25, 7, 7, 1e-6);
[fieldt3u1a, fieldw3u1a, psi3u1a, relE3u1a, conv3u1a, niter3u1a, mallniterc3u1a, Jterms3u1a, maxgrad3u1a, alpha3u1a, invHess3u1a] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1a, options4a, 800, 0.25, 7, 7, 1e-6);
size(conv3u1)
size(conv3u1a)
16-Jterms3u1a.Jmax
16-Jterms3u1.Jmax
figure
plot(0:pi/800:pi, fieldw3u1a(1:801))
hold on
plot(0:pi/800:pi, fieldw3u1(1:801))
figure
plot(0:0.25:800, fieldt3u1)
hold on
plot(0:0.25:800, fieldt3u1a)
[fieldt3u1a, fieldw3u1a, psi3u1a, relE3u1a, conv3u1a, niter3u1a, mallniterc3u1a, Jterms3u1a, maxgrad3u1a, alpha3u1a, invHess3u1a] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(10*(w-0.3))).*sin(2*pi*w/0.3), filterE3u1a, options4a, 800, 0.25, 7, 7, 1e-6);
[fieldt3u1b, fieldw3u1b, psi3u1b, relE3u1b, conv3u1b, niter3u1b, mallniterc3u1b, Jterms3u1b, maxgrad3u1b, alpha3u1b, invHess3u1b] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1b, options4a, 800, 0.25, 7, 7, 1e-6);
filterE3u1b = @(w)20*0.5*(1-tanh(20*(w-0.25)))
options4a.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterE3u1b);
[fieldt3u1b, fieldw3u1b, psi3u1b, relE3u1b, conv3u1b, niter3u1b, mallniterc3u1b, Jterms3u1b, maxgrad3u1b, alpha3u1b, invHess3u1b] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.03*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1b, options4a, 800, 0.25, 7, 7, 1e-6);
16-Jterms3u1a.Jmax
16-Jterms3u1b.Jmax
plot(0:0.25:800, fieldt3u1b)
max(abs(fieldt3u1b-fieldt3u1a))
[fieldt3u1b, fieldw3u1b, psi3u1b, relE3u1b, conv3u1b, niter3u1b, mallniterc3u1b, Jterms3u1b, maxgrad3u1b, alpha3u1b, invHess3u1b] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options4a, 800, 0.25, 7, 7, 1e-6);
16-Jterms3u1b.Jmax
plot(0:0.25:800, fieldt3u1b)
figure
plot(0:0.25:800, fieldt3u1b)
conv3u1b(end)
conv3u1a(end)
max(abs(fieldt3u1b-fieldt3u1b(end:-1:1)))
16-Jterms3u1b.Jmax
16-Jterms3u1a.Jmax
options4b = optionsOCqn(1e-7, 1e4);
options4b.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterE3u1a);
[fieldt3u1b1, fieldw3u1b1, psi3u1b1, relE3u1b1, conv3u1b1, niter3u1b1, mallniterc3u1b1, Jterms3u1b1, maxgrad3u1b1, alpha3u1b1, invHess3u1b1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], fieldw3u1b, filterE3u1a, options4b, 800, 0.25, 7, 7, 1e-8
options4b.invHess0 = invHess3u1b;
[fieldt3u1b1, fieldw3u1b1, psi3u1b1, relE3u1b1, conv3u1b1, niter3u1b1, mallniterc3u1b1, Jterms3u1b1, maxgrad3u1b1, alpha3u1b1, invHess3u1b1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], fieldw3u1b, filterE3u1a, options4b, 800, 0.25, 7, 7, 1e-8);
figure
plot(0:0.25:800, fieldt3u1b1)
max(abs(fieldt3u1b1-fieldt3u1b1(end:-1:1)))
figure
plot(0:pi/800:pi, fieldw3u1b1(1:801))
whos
options4c = options4b;
options4c.invHess0 = invHess3u1a;
[fieldt3u1a1, fieldw3u1a1, psi3u1a1, relE3u1a1, conv3u1a1, niter3u1a1, mallniterc3u1a1, Jterms3u1a1, maxgrad3u1a1, alpha3u1a1, invHess3u1a1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1, 1, fcouplingOp3u1, [-2 6], fieldw3u1a, filterE3u1a, options4c, 800, 0.25, 7, 7, 1e-8);
figure
plot(0:0.25:800, fieldt3u1a1)
whos
clear invHess invHess1 invHess1a invHess2 invHess3u1 invHess3u1a invHess3u1a1 invHess3u1b invHess3u1b1
whos
save HTLSs_harmonic
H3LS_1h = kron(eye(4), H3LS_1);
H3LS_2h = kron(eye(4), H3LS_2);
H3LS_2hs = sparse(kron(eye(4), H3LS_2))
H03
H3LS_2hu = zeros(13);
H3LS_2hu(2:4, 2:4) = H3LS_2h(isingle3, isingle3);
H3LS_2hu(5:7, 5:7) = H3LS_2h(isingle3, isingle3);
H3LS_2hu(8:13, 8:13) = H3LS_2h(idouble3, idouble3)
H3LS_2hus = sparse(H3LS_2hu)
figure
plot(0:0.25:800, psi3u1b1(1,:).*conj(psi3u1b1(1,:)))
plot(0:0.25:800, psi3u1b1(2:4,:).*conj(psi3u1b1(2:4,:)))
plot(0:0.25:800, psi3u1b1(5:7,:).*conj(psi3u1b1(5:7,:)))
plot(0:0.25:800, psi3u1b1(8:13,:).*conj(psi3u1b1(8:13,:)))
max(abs(psi3u1b1(8:13,:).*conj(psi3u1b1(8:13,:))), [], 2)
idouble
idouble3
options5 = optionsOCqn(1e-4, 1e4);
options4.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
[Hoperations3u1c, fcouplingOp3u1c] = Hmats2Hops2(H3harc1us, H3LS_1hus, H3LS_2hus);
[fieldt3u1c, fieldw3u1c, psi3u1c, relE3u1c, conv3u1c, niter3u1c, mallniterc3u1c, Jterms3u1c, maxgrad3u1c, alpha3u1c, invHess3u1c] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6
options5.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
options4.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterE3u1);
[fieldt3u1c, fieldw3u1c, psi3u1c, relE3u1c, conv3u1c, niter3u1c, mallniterc3u1c, Jterms3u1c, maxgrad3u1c, alpha3u1c, invHess3u1c] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6);
size(H3LS_2hus);
size(H3LS_2hus)
[fieldt3u1c, fieldw3u1c, psi3u1c, relE3u1c, conv3u1c, niter3u1c, mallniterc3u1c, Jterms3u1c, maxgrad3u1c, alpha3u1c, invHess3u1c] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6);
16-Jterms3u1c.Jmax
figure
plot(0:0.25:800, fieldt3u1c)
size(fieldw3u1c)
figure
plot(0:0.25:800, psi3u1c(8:13,:).*conj(psi3u1c(8:13,:)))
16-Jterms3u1b1.Jmax
coupling1
max(abs(psi3u1c(8:13,:).*conj(psi3u1c(8:13,:))), [], 2)
800/(10*pi)
plot(0:0.25:800, psi3u1c(2:4,:).*conj(psi3u1c(2:4,:)))
plot(0:0.25:800, psi3u1c(8:13,:).*conj(psi3u1c(8:13,:)))
options6 = optionsOCqn(1e-7, 1e4);
options6.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
options6.invHess0 = invHess3u1c;
[fieldt3u1c1, fieldw3u1c1, psi3u1c1, relE3u1c1, conv3u1c1, niter3u1c1, mallniterc3u1c1, Jterms3u1c1, maxgrad3u1c1, alpha3u1c1, invHess3u1c1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], fieldw3u1c, filterE3u1a, options6, 800, 0.25, 7, 7, 1e-8);
figure
plot(0:0.25:800, fieldt3u1c1)
hold on
plot(0:0.25:800, fieldt3u1c)
[fieldt3u1c1, fieldw3u1c1, psi3u1c1, relE3u1c1, conv3u1c1, niter3u1c1, mallniterc3u1c1, Jterms3u1c1, maxgrad3u1c1, alpha3u1c1, invHess3u1c1] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], [fieldw3u1b1; zeros(1, 3201)], filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6);
16-Jterms3u1c1.Jmax
16-Jterms3u1c.Jmax
conv3u1c(end)
conv3u1c1(end)
conv3u1c(end)-conv3u1c1(end)
figure
plot(0:0.25:800, fieldt3u1c)
hold on
plot(0:0.25:800, fieldt3u1c1)
niter3u1c
niter3u1c1
niter3u1b
366+452
H03
coupling1
target3LSu
%-- 24/09/2020 11:45 --%
load HTLSs_harmonic
H3LS_1h = kron(eye(4), H3LS_1);
H3LS_2h = kron(eye(4), H3LS_2);
H3LS_2hs = sparse(kron(eye(4), H3LS_2))
H03
H3LS_2hu = zeros(13);
H3LS_2hu(2:4, 2:4) = H3LS_2h(isingle3, isingle3);
H3LS_2hu(5:7, 5:7) = H3LS_2h(isingle3, isingle3);
H3LS_2hu(8:13, 8:13) = H3LS_2h(idouble3, idouble3)
H3LS_2hus = sparse(H3LS_2hu)
options5 = optionsOCqn(1e-4, 1e4);
options4.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
[Hoperations3u1c, fcouplingOp3u1c] = Hmats2Hops2(H3harc1us, H3LS_1hus, H3LS_2hus);
load HTLSs_harmonic
options4.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
[Hoperations3u1c, fcouplingOp3u1c] = Hmats2Hops2(H3harc1us, H3LS_1hus, H3LS_2hus);
options5.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
options4.f_max_alpha = get_f_max_alphaOCf(5, 0.25, 800, filterE3u1);
save HTLSs_harmonic
[fieldt3u1c, fieldw3u1c, psi3u1c, relE3u1c, conv3u1c, niter3u1c, mallniterc3u1c, Jterms3u1c, maxgrad3u1c, alpha3u1c, invHess3u1c] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6);
whos
figure
plot(0:0.25:800, fieldt3u1c)
whos
figure
plot(0:pi/800:pi, fieldw3u1c(1:801))
plot(0:pi/800:pi, fieldw3u1c)
plot(0:pi/800:pi, fieldw3u1c(:,1:801))
fieldw3u1c_noise = fieldw3u1c + (rand - 0.5)*max
noise = 0.5*(1-tanh(20*(w3u1-0.25)))*(rand(2, 3201)-0.5)*1e-2*0.02;
size(w3u1)
noise = 0.5*(1-tanh(20*([w3u1;w3u1]-0.25))).*(rand(2, 3201)-0.5)*1e-2*0.02;
figure
plot(w3u1, noise)
fieldw3u1c_noise = fieldw3u1c + noise;
options6 = optionsOCqn(1e-4, 1);
options6.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
[fieldt3u1n, fieldw3u1n, psi3u1n, relE3u1n, conv3u1n, niter3u1n, mallniterc3u1n, Jterms3u1n, maxgrad3u1n, alpha3u1n, invHess3u1n] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], fieldw3u1c_noise, filterE3u1a, options6, 800, 0.25, 7, 7, 1e-6);
options6 = optionsOCqn(1e-4, 0);
options6.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 800, filterE3u1a, 2);
[fieldt3u1n, fieldw3u1n, psi3u1n, relE3u1n, conv3u1n, niter3u1n, mallniterc3u1n, Jterms3u1n, maxgrad3u1n, alpha3u1n, invHess3u1n] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], fieldw3u1c_noise, filterE3u1a, options6, 800, 0.25, 7, 7, 1e-6);
Jterms3u1n
16-Jterms3u1n.Jmax
max(max(fieldwu31c))
max(max(fieldw3u1c))
(max(fieldw3u1c), [], 2)
max(fieldw3u1c, [], 2)
noise = 0.5*(1-tanh(20*([w3u1;w3u1]-0.25))).*(rand(2, 3201)-0.5)*1e-2.*max(abs(fieldw3u1c, [], 2));
noise = 0.5*(1-tanh(20*([w3u1;w3u1]-0.25))).*(rand(2, 3201)-0.5)*1e-2.*max(abs(fieldw3u1c), [], 2);
fieldw3u1c_noise = fieldw3u1c + noise;
[fieldt3u1n, fieldw3u1n, psi3u1n, relE3u1n, conv3u1n, niter3u1n, mallniterc3u1n, Jterms3u1n, maxgrad3u1n, alpha3u1n, invHess3u1n] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1c, 2, fcouplingOp3u1c, [-2 6], fieldw3u1c_noise, filterE3u1a, options6, 800, 0.25, 7, 7, 1e-6);
16-Jterms3u1n.Jmax
16-Jterms3u1c.Jmax
(16-Jterms3u1c.Jmax)/16
(16-Jterms3u1n.Jmax)/16
whos
clear invHess3u1c invHess3u1n
save HTLSs_harmonic
%-- 14/10/2020 12:08 --%
load HTLSs_harmonic
H3LS_1
H3LS_2
H3LS_1 = kron(H03e, eye(3))
H03
H03e = diag(0:2)
H3LSe_1 = kron(H03e, eye(3))
H3LSe_2 = kron(eye(3), H03e)
H3LSe_1h = kron(eye(4), H3LSe_1);
H3LSe_2h = kron(eye(4), H3LSe_2);
H3LS_1hu = zeros(13);
H3LSe_1hu = zeros(13);
H3LS_1hu(2:4, 2:4) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(5:7, 5:7) = H3LS_1h(isingle3, isingle3);
H3LS_1hu(8:13, 8:13) = H3LS_1h(idouble3, idouble3)
H3LSe_1hu(2:4, 2:4) = H3LSe_1h(isingle3, isingle3);
H3LSe_1hu(5:7, 5:7) = H3LSe_1h(isingle3, isingle3);
H3LSe_1hu(8:13, 8:13) = H3LSe_1h(idouble3, idouble3)
H3LSe_2 = kron(eye(3), H03)
H3LSe_2 = kron(eye(3), H03e)
H3LSe_2h = kron(eye(4), H3LSe_2);
H3LSe_2hu = zeros(13);
H3LSe_2hu(2:4, 2:4) = H3LSe_2h(isingle3, isingle3);
H3LSe_2hu(5:7, 5:7) = H3LSe_2h(isingle3, isingle3);
H3LSe_2hu(8:13, 8:13) = H3LSe_2h(idouble3, idouble3)
H3LSe_2hus = sparse(H3LSe_2hu)
[fieldt3u1ce, fieldw3u1ce, psi3u1ce, relE3u1ce, conv3u1ce, niter3u1ce, mallniterc3u1ce, Jterms3u1ce, maxgrad3u1ce, alpha3u1ce, invHess3u1ce] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6);
[Hoperations3u1ce, fcouplingOp3u1ce] = Hmats2Hops2(H3harc1us, H3LSe_1hus, H3LSe_2hus);
H3harc1us
H3harcus
H3harc1us = sparse(H3harc1u)
H3harc1
load HTLSs_harmonic
H3harc1
load HTLSs_harmonic
H3harc1
[Hoperations3u1ce, fcouplingOp3u1ce] = Hmats2Hops2(H3harc1us, H3LSe_1hus, H3LSe_2hus);
H3LSe_1hus = sparse(H3LSe_1hu)
H3LSe_2hus = sparse(H3LSe_2hu)
[Hoperations3u1ce, fcouplingOp3u1ce] = Hmats2Hops2(H3harc1us, H3LSe_1hus, H3LSe_2hus);
[fieldt3u1ce, fieldw3u1ce, psi3u1ce, relE3u1ce, conv3u1ce, niter3u1ce, mallniterc3u1ce, Jterms3u1ce, maxgrad3u1ce, alpha3u1ce, invHess3u1ce] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 800, 0.25, 7, 7, 1e-6);
16-Jterms3u1ce.Jmax
16-Jterms3u1c.Jmax
figure
plot(0:0.25:800, fieldt3u1c)
figure
plot(0:0.25:800, fieldt3u1ce)
hold on
plot(0:0.25:800, fieldt3u1c)
figure
plot(0:0.25:800, fieldt3u1ce-fieldt3u1c)
figure
plot((0:0.25:800)/(5*2*pi), fieldt3u1ce*5)
plot(0:0.25:800, fieldt3u1ce)
xlabel('$\omega_0t$', 'interpreter', 'latex')
ylabel('$\eta_n(\omega_0t)/\omega_0$',  'interpreter', 'latex')
800/(10*pi)
clear invHess3u1ce
save HTLSs_harmonic
%-- 18/10/2020 14:23 --%
load HTLSs_harmonic
800/(10*pi)
10*pi/3
tdiscrete = 10*pi/6:10*pi/3:800;
size(tdiscrete)
w
size(fieldw3u1ce)
t800 = 0:pi/800:pi/0.25;
t800 = 0:0.25:800;
w800 = 0:pi/800:pi/0.25;
dctfactor = 800/(sqrt(3200*pi));
fieldt3u1ce_dis = sqrt(2/pi)*pi/800*(cos(tdiscrete.'*w800)*fieldw3u1ce.').';
size(fieldt3u1ce_dis)
figure
plot(0:0.25:800, fieldt3u1ce)
hold on
plot(tdiscrete, fieldt3u1ce_dis)
figure
plot(tdiscrete, fieldt3u1ce_dis(1,:))
hold on
plot(0:0.25:800, fieldt3u1ce(1, :))
fieldt3u1ce_dis = sqrt(2/pi)*pi/800*(cos(tdiscrete.'*w800)*[0.5*fieldw3u1ce(1), fieldw3u1ce(2:3200), 0.5*fieldw3u1ce(3201)].').';
figure
plot(0:0.25:800, fieldt3u1ce)
hold on
plot(tdiscrete, fieldt3u1ce_dis)
fieldt3u1ce_dis = sqrt(2/pi)*pi/800*(cos(tdiscrete.'*w800)*[0.5*fieldw3u1ce(1); fieldw3u1ce(2:3200).'; 0.5*fieldw3u1ce(3201)]).';
plot(tdiscrete, fieldt3u1ce_dis)
fieldt3u1ce_dis = sqrt(2/pi)*pi/800*(cos(tdiscrete.'*w800)*[0.5*fieldw3u1ce(:,1), fieldw3u1ce(:, 2:3200), 0.5*fieldw3u1ce(:, 3201)].').';
plot(tdiscrete, fieldt3u1ce_dis)
save HTLSs_harmonic
%-- 19/10/2020 15:29 --%
load HTLSs_harmonic
10*pi/3*76
tdiscrete(end)
tdiscrete(end)+10*pi/6
testPWCharmonic
Hoperations3u1ce
Hop3u1ce = @(v, field) H3harc1us - field(1)*(H3LSe_1hus*v)-field(2)*(H3LSe_2hus*v)
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, tdiscrete, [-2, 6], 10*pi/3*76, 76, 7);
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce, [-2, 6], 10*pi/3*76, 76, 7);
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 7);
size(H3harc1us), size(H3LSe_1hus), size(H3LSe_2hus)
Hop3u1ce = @(v, field) H3harc1us*v - field(1)*(H3LSe_1hus*v)-field(2)*(H3LSe_2hus*v)
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 7);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 15);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 30);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 50);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 100);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 200);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce_dis, [-2, 6], 10*pi/3*76, 76, 100);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
[psi3u1ce(:,end).*conj(psi3u1ce(:,end)), psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))]
figure
plot(tdiscrete, fieldt3u1ce_dis)
figure
plot(0:0.25:800, conj(psi3u1ce).*psi3u1ce)
plot(0:0.25:800, conj(psi3u1ce(2:4, :)).*psi3u1ce(2:4, :))
hold on
plot(tdiscrete, conj(psiPWC3u1ce(2:4, :)).*psiPWC3u1ce(2:4, :))
size(psiPWC3u1ce)
plot(0:10*pi/3:76*10*pi/3, conj(psiPWC3u1ce(2:4, :)).*psiPWC3u1ce(2:4, :))
10*pi/3
10*pi/6
3200/40
tdiscrete1 = 5:10:800;
size(tdiscrete(1))
size(tdiscrete1)
tdiscrete1(end)
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce(:, 1:20:end), [-2, 6], 800, 80, 100);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
size(psiPWC3u1ce)
psiPWC3u1ce = SchrPWCcheb(Hop3u1ce, u0_3, fieldt3u1ce(:, 21:40:3201), [-2, 6], 800, 80, 100);
psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))
[psi3u1ce(:,end).*conj(psi3u1ce(:,end)), psiPWC3u1ce(:,end).*conj(psiPWC3u1ce(:,end))]
psiPWC3u1ce(:,end)'.*ftarget3u1(psiPWC3u1ce(:,end))
psiPWC3u1ce(:,end)'*ftarget3u1(psiPWC3u1ce(:,end))
16-psiPWC3u1ce(:,end)'*ftarget3u1(psiPWC3u1ce(:,end))
(16-psiPWC3u1ce(:,end)'*ftarget3u1(psiPWC3u1ce(:,end)))/16
16-Jterms3u1ce.Jmax
(16-Jterms3u1ce.Jmax)/16
save HTLSs_harmonic
%-- 21/10/2020 10:04 --%
load HTLSs_harmonic
H0harc
Hharmonic
Hharmonic2 = Hharmonic(1:3,1:3)*1.3
H02modes = kron(Hharmonic2, eye(27)) + kron(eye(3), H3harc(1:27,1:27));
size(H02modes)
H02modess = sparse(kron(Hharmonic2, eye(27)) + kron(eye(3), H3harc(1:27,1:27)))
diag(H02modess)
Hmode2couplings = sparse(coupling*(kron(adag3, kron(eye(3), sigmam31)) + kron(a3, kron(eye(3), sigmap31)) - kron(adag3, kron(eye(3), sigmam32)) - kron(a3, kron(eye(3), sigmap32))))
Hmode2couplings = sparse(coupling1*(kron(adag3, kron(eye(3), sigmam31)) + kron(a3, kron(eye(3), sigmap31)) - kron(adag3, kron(eye(3), sigmam32)) - kron(a3, kron(eye(3), sigmap32))))
H02modess = sparse(kron(Hharmonic2, eye(27)) + kron(eye(3), H3harc1(1:27,1:27)))
H2modess = H02modess + Hmode2coupling;
H2modess = H02modess + Hmode2couplings;
H2modess = H02modess + Hmode2couplings
Htls
H0TLSs
HTLSs
HTLSs1 = H0TLSs + coupling1*Hcoupling
UsymTLSs
psiTLS1 = fMdiag(HTLSs1, @(H) exp(-1i*H), psi0, 2*pi*30, 1e3);
psi0
figure
plot(0:(2*pi*30/1e3):2*pi*30)
plot(0:(2*pi*30/1e3):2*pi*30, psiTLS1.*conj(psiTLS1))
pi/4/(0.06)
pi/0.06
pi/2/(0.03)
800/ans
pi/2/(0.03)
pi/2/(0.03)/(10*pi)
[fieldt3u2, fieldw3u2, psi3u2, relE3u2, conv3u2, niter3u2, mallniterc3u2, Jterms3u2, maxgrad3u2, alpha3u2, invHess3u2] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u1a, options5, 200, 0.25, 7, 7, 1e-6
options7 = optionsOCqn(1e-4, 1e4);
options7.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 200, filterE3u1a, 2
size(filterE3u1a)
filterE3u1a = @(w)10*0.5*(1-tanh(20*(w-0.25)))
filterE3u2 = @(w)2.5*0.5*(1-tanh(20*(w-0.25)))
size(psi3u2)
options7.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 200, filterE3u2, 2);
[fieldt3u2, fieldw3u2, psi3u2, relE3u2, conv3u2, niter3u2, mallniterc3u2, Jterms3u2, maxgrad3u2, alpha3u2, invHess3u2] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u2, options7, 200, 0.25, 7, 7, 1e-6);
16-Jterms3u2.Jmax
figure
plot(0:0.25:200, fieldt3u2)
size(fieldt3u2)
max(abs(fieldtu2(1,:)-fieldtu2(2,:)))
max(abs(fieldt3u2(1,:)-fieldt3u2(2,:)))
Hoperations3u1ce
H3LSe_1hus
H3LSe_2hus
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), filterE3u3, options7, 400, 0.25, 7, 7, 1e-6
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options7, 400, 0.25, 7, 7, 1e-6);
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options5, 400, 0.25, 7, 7, 1e-6);
options8 = optionsOCqn(1e-4, 1e4);
options3u3 = optionsOCqn(1e-4, 1e4);
options3u3.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 200, @(w)5*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
options3u3.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 400, @(w)5*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
options3u3.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 400, @(w)5*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
16-Jterms3u3.Jmax
figure
plot(0:0.25:400, fieldt3u3)
max(abs(fieldt3u3-fieldt3u3(:,end:-1:1)),[],2)
figure
plot(0:0.25:800, fieldt3u1ce)
figure
plot(0:0.25:800, sum(fieldt3u1ce))
hold on
plot(0:0.25:800, sum(fieldt3u1ce))
plot(0:0.25:800, sum(fieldt3u1ce(end:-1:1)))
clf
plot(0:0.25:800, sum(fieldt3u1ce(:,end:-1:1)))
clf
plot(0:0.25:800, sum(fieldt3u1ce))
hold on
plot(0:0.25:800, sum(fieldt3u1ce(:,end:-1:1)))
options3u4.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 300, @(w)3.75*0.5*(1-tanh(20*(w-0.25))), 2);
options3u4 = optionsOCqn(1e-4, 1e4);
options3u4.f_max_alpha = get_f_max_alphaOCf_multiE(5, 0.25, 300, @(w)3.75*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u4, fieldw3u4, psi3u4, relE3u4, conv3u4, niter3u4, mallniterc3u4, Jterms3u4, maxgrad3u4, alpha3u4, invHess3u4] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)3.75*0.5*(1-tanh(20*(w-0.25))), options3u4, 300, 0.25, 7, 7, 1e-6);
16-Jterms3u4.Jmax
figure
16-Jterms3u2.Jmax
(16-Jterms3u2.Jmax)/16
(16-Jterms3u4.Jmax)/16
(16-Jterms3u3.Jmax)/16
plot(0:0.25:300, fieldt3u4)
max(abs(fieldt3u4-fieldt3u4(:,end:-1:1)),[],2)
hold on
plot(0:0.25:300, fieldt3u4(:, end:-1:1))
mean(fieldt3u4,[],2)
mean(fieldt3u4,2)
mean(fieldt3u2,2)
whos
clear invHess3u2 invHess3u3 invHess3u4
save HTLSs_harmonic
%-- 08/11/2020 12:08 --%
load HTLSs_harmonic
size(target3LSu(:))
size(target3LSu)
target3LSu
whos
H3harc1us
H3harc1usI = H3harc1us - diag([0, ones(6, 1), 2*ones(6,1)])
H3harc1usI = H3harc1us - diag([0; ones(6, 1); 2*ones(6,1)])
H3harc1usI = H3harc1us - sparse(diag([0; ones(6, 1); 2*ones(6,1)]))
nnz(H3harc1us)
nnz(H3harc1usI)
diag(H3harc1usI)
eig(H3harc1usI)
H3LSe_1hus
H3LSe_2hus
eig(H3harc1usI + 0.05*(H3LSe_1hus + H3LSe_2hus))
eig(H3harc1usI + 0.1*(H3LSe_1hus + H3LSe_2hus))
eig(H3harc1usI - 0.1*(H3LSe_1hus + H3LSe_2hus))
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
[Hoperations3u1ceI, fcouplingOp3u1ceI] = Hmats2Hops2(H3harc1usI, H3LSe_1hus, H3LSe_2hus);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
(16-Jterms3u3.Jmax)/16
(16-Jterms3u3I.Jmax)/16
figure
plot(0:0.25:400, fieldt3u3)
plot(0:0.25:400, fieldt3u3I)
plot(0:0.25:400, fieldt3u3)
options3u3I = options3u3
options3u3I.maxNiter = 0;
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 0.25, 7, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
ftarget3u3 = phi2ftarget(exp(1i*)target3LSu(:))
H0I = sparse(diag([0; ones(6, 1); 2*ones(6,1)]))
ftarget3u3 = phi2ftarget(exp(1i*[0; ones(6, 1); 2*ones(6,1)]*400).*target3LSu(:))
clear ftarget3u3
ftarget3u3I = phi2ftarget(exp(1i*[0; ones(6, 1); 2*ones(6,1)]*400).*target3LSu(:))
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 0.25, 7, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
(16-Jterms3u3.Jmax)/16
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 0.25, 7, 7, 1e-6);
(16-Jterms3u3.Jmax)/16
(16-Jterms3u3I.Jmax)/16
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 5, 7, 7, 1e-6);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 4, 7, 7, 1e-6);
figure
plot(0:0.25:400, fieldt3u3)
hold on
plot(0:4:400, fieldt3u3I)
max(abs(fieldt3u3(1:16:end)-fieldt3u3I))
max(abs(fieldt3u3(:, 1:16:end)-fieldt3u3I), [], 2)
mallniterc3u3I
mallniterc3u3
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 4, 7, 7, 1e-6);
[fieldt3u3, fieldw3u3, psi3u3, relE3u3, conv3u3, niter3u3, mallniterc3u3, Jterms3u3, maxgrad3u3, alpha3u3, invHess3u3] = OClimf_qn(u0_3, ftarget3u1, Hoperations3u1ce, 2, fcouplingOp3u1ce, [-2 6], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 0.25, 7, 7, 1e-6);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 4, 7, 7, 1e-6);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 7, 7, 1e-6);
mallniterc3u3I
max(abs(fieldt3u3(:, 1:10:end)-fieldt3u3I), [], 2)
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 7, 7, 1e-6);
146*pi/400
pi/2.5
pi/4
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 5, 7, 7, 1e-6);
max(abs(fieldt3u3(:, 1:20:end)-fieldt3u3I), [], 2)
mallniterc3u3I
mallniterc3u3
6/400*2*pi
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], fieldw3u3, @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 2.5, 7, 7, 1e-6);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 2.5, 7, 7, 1e-6);
options3u3I = options3u3;
options3u3I.f_max_alpha = get_f_max_alphaOCf_multiE(1, 2.5, 400, @(w)5*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3, 400, 2.5, 7, 7, 1e-6);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 7, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
figure
plot(0:2.5:400, fieldt3u3I)
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 9, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
figure
plot(0:2.5:400, fieldt3u3I)
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 1, 9, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
figure
plot(0:2.5:400, fieldt3u3I)
plot(0:1:400, fieldt3u3I)
options3u3I.f_max_alpha = get_f_max_alphaOCf_multiE(0.2, 2.5, 400, @(w)5*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 1, 9, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
(16-Jterms3u3.Jmax)/16
figure
plot(0:1:400, fieldt3u3I)
max(abs(fieldt3u3(:, 1:4:end)-fieldt3u3I), [], 2)
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.2 0.2], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 9, 7, 1e-6);
(16-Jterms3u3.Jmax)/16
(16-Jterms3u3I.Jmax)/16
figure
plot(0:2.5:400, fieldt3u3I)
max(abs(fieldt3u3(:, 1:10:end)-fieldt3u3I), [], 2)
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.3 0.3], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 9, 7, 1e-6);
%-- 10/11/2020 10:31 --%
options3u3I.f_max_alpha = get_f_max_alphaOCf_multiE(0.2, 2.5, 400, @(w)5*0.5*(1-tanh(20*(w-0.25))), 2);
[fieldt3u3I, fieldw3u3I, psi3u3I, relE3u3I, conv3u3I, niter3u3I, mallniterc3u3I, Jterms3u3I, maxgrad3u3I, alpha3u3I, invHess3u3I] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1ceI, 2, fcouplingOp3u1ce, [-0.3 0.3], @(w)0.05*0.5*(1-tanh(20*(w-0.25))).*sin(2*pi*w/0.25), @(w)5*0.5*(1-tanh(20*(w-0.25))), options3u3I, 400, 2.5, 9, 7, 1e-6);
(16-Jterms3u3I.Jmax)/16
max(abs(fieldt3u3([2;1], 1:16:end)-fieldt3u3I), [], 2)
max(abs(fieldt3u3(:, 1:10:end)-fieldt3u3I), [], 2)
whos
clear mx U error invHess3u3I
whos
clear t tmiddle
save HTLSs_harmonic
H3harc1usIg = H3harc1usI/0.03
Hoperations3u1cIg = Hmats2Hops2(H3harc1usIg, H3LSe_1hus, H3LSe_2hus);
options3u3Ig = optionsOCqn(1e-4, 1e4);
options3u3Ig.f_max_alpha = get_f_max_alphaOCf_multiE(6, 0.075, 12, @(w)5*0.5*(1-tanh(20*(w-0.25/0.03))), 2);
options3u3Ig.maxNiter = 0;
[fieldt3u3Ig, fieldw3u3Ig, psi3u3Ig, relE3u3Ig, conv3u3Ig, niter3u3Ig, mallniterc3u3Ig, Jterms3u3Ig, maxgrad3u3Ig, alpha3u3Ig] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], fieldw3u3I, @(w)5*0.5*(1-tanh(0.6*(w-0.25/0.03))), options3u3Ig, 12, 0.075, 9, 7, 1e-6);
max(max(abs(psi3u3I-psi3u3Ig)))
save HTLSs_harmonic
options3u3Ig.maxNiter = 1e4;
[fieldt3u3Ig, fieldw3u3Ig, psi3u3Ig, relE3u3Ig, conv3u3Ig, niter3u3Ig, mallniterc3u3Ig, Jterms3u3Ig, maxgrad3u3Ig, alpha3u3Ig] = OClimf_qn(u0_3, ftarget3u3I, Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)5/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options3u3Ig, 12, 0.075, 9, 7, 1e-6);
(16-Jterms3u3Ig.Jmax)/16
max(abs(fieldt3u3Ig-fieldt3u3I/0.03), [], 2)
5/0.03
options_gate = optionsOCqn(1e-4, 1e4);
options_gate.f_max_alpha = get_f_max_alphaOCf_multiE(0.2/0.03, 0.075, 12, @(w)5/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
[fieldtg, fieldwg, psig, relEg, convg, niterg, mallnitercg, Jtermsg, maxgradg, alphag] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate, 12, 0.075, 9, 7, 1e-6);
[fieldtg, fieldwg, psig, relEg, convg, niterg, mallnitercg, Jtermsg, maxgradg, alphag] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)40/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate, 12, 0.075, 9, 7, 1e-6);
[fieldtg, fieldwg, psig, relEg, convg, niterg, mallnitercg, Jtermsg, maxgradg, alphag] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate, 12, 0.075, 9, 7, 1e-6);
Jtermsg
[fieldtg1, fieldwg1, psig1, relEg1, convg1, niterg1, mallnitercg1, Jtermsg1, maxgradg1, alphag1] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], fieldwg, @(w)40/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate, 12, 0.075, 9, 7, 1e-6);
figure
plot(0:0.075:12, fieldtg)
plot(0:0.075:12, fieldtg1)
Jtermsg1
[fieldtg2, fieldwg2, psig2, relEg2, convg2, niterg2, mallnitercg2, Jtermsg2, maxgradg2, alphag2] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], fieldw3u3Ig, @(w)80/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate, 12, 0.075, 9, 7, 1e-6);
1-Jtermsg2.Jmax
1-Jtermsg1.Jmax
figure
plot(0:0.075:12, fieldtg2)
hold on
plot(0:0.075:12, fieldt3u3Ig)
Jtermsg2
phase2 = (target3LSu(2:4)'*psig1(2:4,end))/abs((target3LSu(2:4)'*psig1(2:4,end)))
phase3 = (target3LSu(5:7)'*psig1(5:7,end))/abs((target3LSu(5:7)'*psig1(5:7,end)))
phase4 = (target3LSu(8:13)'*psig1(8:13,end))/abs((target3LSu(8:13)'*psig1(8:13,end)))
1 - conj(phase2*phase3)*phase4
conj(phase2*phase3)*phase4
real((target3LSu(8:13)'*psig1(8:13,end))*conj(phase2*phase3))^2
(target3LSu(2:4)'*psig(2:4,end))*conj(target3LSu(2:4)'*psig(2:4,end))
(target3LSu(5:7)'*psig(5:7,end))*conj(target3LSu(5:7)'*psig(5:7,end))
(target3LSu(8:13)'*psig(8:13,end))*conj(target3LSu(8:13)'*psig(8:13,end))
(target3LSu(2:4)'*psig1(2:4,end))*conj(target3LSu(2:4)'*psig1(2:4,end))
(target3LSu(5:7)'*psig1(5:7,end))*conj(target3LSu(5:7)'*psig1(5:7,end))
(target3LSu(8:13)'*psig1(8:13,end))*conj(target3LSu(8:13)'*psig1(8:13,end))
Jtermsg2
Jterms3u3Ig
Jterms3u3Ig/16
Jterms3u3Ig.Jenergy/16
(1 + (target3LSu(2:4)'*psig1(2:4,end))*conj(target3LSu(2:4)'*psig1(2:4,end))
+ real((target3LSu(8:13)'*psig1(8:13,end))*conj(phase2*phase3))^2)
(1  + real((target3LSu(8:13)'*psig1(8:13,end))*conj(phase2*phase3))^2)
1 - (1 + (target3LSu(2:4)'*psig1(2:4,end))*conj(target3LSu(2:4)'*psig1(2:4,end))
+ real((target3LSu(8:13)'*psig1(8:13,end))*conj(phase2*phase3))^2)
1 - (1 + (target3LSu(2:4)'*psig1(2:4,end))*conj(target3LSu(2:4)'*psig1(2:4,end) + (target3LSu(5:7)'*psig1(5:7,end))*conj(target3LSu(5:7)'*psig1(5:7,end)) + real((target3LSu(8:13)'*psig1(8:13,end))*conj(phase2*phase3))^2)/4
1 - (1 + (target3LSu(2:4)'*psig1(2:4,end))*conj(target3LSu(2:4)'*psig1(2:4,end)) + (target3LSu(5:7)'*psig1(5:7,end))*conj(target3LSu(5:7)'*psig1(5:7,end)) + real((target3LSu(8:13)'*psig1(8:13,end))*conj(phase2*phase3))^2)/4
1-Jtermsg1.Jmax
[Jmax_overlap, phases] = Uoverlap_gate(psi(:, end), target3LSu(:), [1 4 7])
[Jmax_overlap, phases] = Uoverlap_gate(psig1(:, end), target3LSu(:), [1 4 7])
1- Jmax_overlap
[Jmax_overlapg1, phasesg1] = Uoverlap_gate(psig1(:, end), target3LSu(:), [1 4 7])
[Jmax_overlapg2, phasesg2] = Uoverlap_gate(psig2(:, end), target3LSu(:), [1 4 7])
(16-Jterms3u3Ig.Jmax)/16
1- Jmax_overlapg2
1 -
1-Jtermsg1.Jmax
1-Jtermsg2.Jmax
2*pi/0.075
6/0.075
[fieldtg3, fieldwg3, psig3, relEg3, convg3, niterg3, mallnitercg3, Jtermsg3, maxgradg3, alphag3, invHessg3] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
options_gate3 = optionsOCqn(1e-4, 1e4);
options_gate3.f_max_alpha = get_f_max_alphaOCf_multiE(0.2/0.03, 0.075, 6, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
[fieldtg3, fieldwg3, psig3, relEg3, convg3, niterg3, mallnitercg3, Jtermsg3, maxgradg3, alphag3, invHessg3] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
[fieldtg3, fieldwg3, psig3, relEg3, convg3, niterg3, mallnitercg3, Jtermsg3, maxgradg3, alphag3, invHessg3] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)80/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
[fieldtg3, fieldwg3, psig3, relEg3, convg3, niterg3, mallnitercg3, Jtermsg3, maxgradg3, alphag3, invHessg3] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)40/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
[fieldtg3, fieldwg3, psig3, relEg3, convg3, niterg3, mallnitercg3, Jtermsg3, maxgradg3, alphag3, invHessg3] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
[fieldtg3a, fieldwg3a, psig3a, relEg3a, convg3a, niterg3a, mallnitercg3a, Jtermsg3a, maxgradg3a, alphag3a, invHessg3a] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], fieldwg3a, @(w)40/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
[fieldtg3a, fieldwg3a, psig3a, relEg3a, convg3a, niterg3a, mallnitercg3a, Jtermsg3a, maxgradg3a, alphag3a, invHessg3a] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations3u1cIg, 2, fcouplingOp3u1ce, [-10 10], fieldwg3, @(w)40/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
figure
plot(0:0.075:6, fieldtg3a)
plot(0:0.075:6, fieldtg3)
plot(0:0.075:6, fieldtg3a)
Jtermsg3
Jtermsg3a
H3harc1usIg
H3harcus10nl = H3harc1usIg;
H3harcus10nl(8, 8) = -10;
H3harcus10nl(10, 10) = -10
eig(H3harcus10nl)
Hoperations10nl = Hmats2Hops2(H3harcus10nl, H3LSe_1hus, H3LSe_2hus);
[fieldtg4, fieldwg4, psig4, relEg4, convg4, niterg4, mallnitercg4, Jtermsg4, maxgradg4, alphag4, invHessg4] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations10nl, 2, fcouplingOp3u1ce, [-15 15], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
figure
plot(0:0.075:6, fieldtg4)
Jtermsg4
[Jmax_overlapg4, phasesg4] = Uoverlap_gate(psig4(:, end), target3LSu(:), [1 4 7])
[fieldtg4a, fieldwg4a, psig4a, relEg4a, convg4a, niterg4a, mallnitercg4a, Jtermsg4a, maxgradg4a, alphag4a, invHessg4a] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations10nl, 2, fcouplingOp3u1ce, [-15 15], fieldwg4, @(w)20/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
figure
plot(0:0.075:6, fieldtg4a)
Jtermsg4a
[Jmax_overlapg4a, phasesg4a] = Uoverlap_gate(psig4(:, end), target3LSu(:), [1 4 7])
[Jmax_overlapg4a, phasesg4a] = Uoverlap_gate(psig4a(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlapg4a
[fieldtg4b, fieldwg4b, psig4b, relEg4b, convg4b, niterg4b, mallnitercg4b, Jtermsg4b, maxgradg4b, alphag4b, invHessg4b] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations10nl, 2, fcouplingOp3u1ce, [-15 15], fieldwg4, @(w)40/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
Jtermsg4b
figure
plot(0:0.075:6, fieldtg4b)
plot(0:0.075:6, fieldtg4)
hold on
plot(0:0.075:6, fieldtg4b)
[fieldtg4c, fieldwg4c, psig4c, relEg4c, convg4c, niterg4c, mallnitercg4c, Jtermsg4c, maxgradg4c, alphag4c, invHessg4c] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations10nl, 2, fcouplingOp3u1ce, [-15 15], fieldwg4, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
Jtermsg4c
1-Jtermsg4c.Jmax
figure
plot(0:0.075:6, fieldtg4c)
[fieldtg4c, fieldwg4c, psig4c, relEg4c, convg4c, niterg4c, mallnitercg4c, Jtermsg4c, maxgradg4c, alphag4c, invHessg4c] = OClimf_gate(u0_3, target3LSu(:), [1 4 7], Hoperations10nl, 2, fcouplingOp3u1ce, [-15 15], fieldwg4, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
1/0.02
[Jmax_overlapg4c, phasesg4c] = Uoverlap_gate(psig4a(:, end), target3LSu(:), [1 4 7])
[Jmax_overlapg4c, phasesg4c] = Uoverlap_gate(psig4c(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlapg4c
figure
plot((0:0.075:6)/(2*pi), fieldtg4c)
xlabel('$gt$', 'interpreter', 'latex')
ylabel('$f_n(gt)/g$',  'interpreter', 'latex')
whos
save HTLSs_harmonic
%-- 13/11/2020 14:01 --%
doc kronecker
C = {[1 2 3; 4 5 6; 7 8 9], [5 6; 7 8], [1;2]}
length(C)
C{1}
C{:}
C{2:3}
cell{C{2:3}}
C(2:3)
Mresult = multi_kron(C)
M = kron(C{1}, kron(C{2}, C{3}))
Mresult -M
Mresult - M
M = kron(C{1}, multi_kron(C{2:3}))
M = kron(C{1}, multi_kron(C(2:3)))
Mresult - M
doc prod
v = [1 2 3 4]
cumprod(v)
cumprod(v, 'reverse')
kron_index = Mkron_index(4, [2 3])
kron_index = Mkron_index([2 4], [4 2 3])
kron_index = kron_index([2 4], [4 2 3])
clear all
kron_index = kron_index([2 4], [4 2 3])
kindex = kron_index([2 4], [4 2 3])
kindex = kron_index([2 4], [4 2 3; 2 1 2])
doc sum
kindex = kron_index([2 4], [4 2 3; 2 1 2])
%-- 16/11/2020 11:42 --%
50/2.4
50*2.4
Hoperations10nl
load HTLSs_harmonic
Hoperations10nl
Hop10nl = @(v, field) H3harcus10nl*v - field(1)*(H3LSe_1hus*v) - field(2)*(H3LSe_2hus*v)
psi4cPWC40 = SchrPWCcheb(Hop10nl, u0_3, fieldtg4c(:, 2:2:end-1), [-15 15], 6, 40, 20);
[Jmax_overlap4cPWC40, phases4cPWC40] = Uoverlap_gate(psi4cPWC40(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlap4cPWC40
psi4cPWC20 = SchrPWCcheb(Hop10nl, u0_3, fieldtg4c(:, 3:4:end-1), [-15 15], 6, 20, 40);
psi4cPWC20(:,end).*conj(psi4cPWC20(:,end))
4 - sqnorm(psi4cPWC20(:,end))
[Jmax_overlap4cPWC20, phases4cPWC20] = Uoverlap_gate(psi4cPWC20(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlap4cPWC20
plot(0:0.075:6, fieldtg4c)
6/20
figure
plot(0:0.075:6, conj(psi4cPWC20).*psi4cPWC20)
plot(2*0.075:4*0.075:6, conj(psi4cPWC20).*psi4cPWC20)
size(2*0.075:4*0.075:6)
plot(0:4*0.075:6, conj(psi4cPWC20).*psi4cPWC20)
plot(0:4*0.075:6, conj(psi4cPWC20(:,2:4)).*psi4cPWC20(:,2:4))
plot(0:4*0.075:6, conj(psi4cPWC20(2:4,:)).*psi4cPWC20(2:4,:))
%-- 17/11/2020 11:07 --%
load HTLSs_harmonic