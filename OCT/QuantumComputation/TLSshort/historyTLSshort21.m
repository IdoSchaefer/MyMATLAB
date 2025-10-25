%-- 14/03/2021 11:44 --%
H0 = [0 0; 0 1]
mu = [0 1; 1 0]
eig(Ho+5*mu)
eig(H0+5*mu)
eig(H0-5*mu)
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) sin(t), 1e-3, options, 2*pi, 2*pi/1e3, 7, 7, 1e-6, 1e4);
options = optionsOCqn(1e-4,1e4);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) sin(t), 1e-3, options, 2*pi, 2*pi/1e3, 7, 7, 1e-4, 1e4);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) sin(t), 1e-3, options, 2*pi, 2*2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) sin(t), 1e-3, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
J1
1-J1
figure
plot(0:2*pi/1e2:2*pi, field)
save TLS2pi
load TLS2pi
imaxfield
imax
whos
imaxf
(2*pi-5.78)
(2*pi-5.78)/pi
imaxf
field_cos = [field(93:end), field(2:93)];
figure
plot(0:2*pi/1e2:2*pi, field_cos)
hold on
plot(0:2*pi/1e2:2*pi, maxfield*cos(0:2*pi/1e2:2*pi))
figure
plot(0:2*pi/1e2:2*pi, field_cos - maxfield*cos(0:2*pi/1e2:2*pi))
delta_cos = field_cos - maxfield*cos(0:2*pi/1e2:2*pi);
figure
plot(0:2*pi/1e2:2*pi, delta_cos + delta_cos(end:-1:1))
figure
plot(0:2*pi/1e2:pi, delta_cos(1:51) + delta_cos(end:-1:51))
plot(0:2*pi/1e2:pi, delta_cos(1:51) + delta_cos(51:101))
[allfield4, field4, psi4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, alpha4, invHess4] = OCqn([1;0], [0;1], H0, [-6 6], mu, field, 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield4, field4, psi4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, alpha4, invHess4] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield, 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield4, field4, psi4, relE4, conv4, niter4, mallniterc4, J14, maxgrad4, alpha4, invHess4] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield.', 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
figure
plot(0:2*pi/1e2:2*pi, field)
hold on
plot(0:2*pi/1e2:2*pi, field4)
max(abs(field-field4))
figure
field4_cos = [field4(93:end), field4(2:93)];
delta_cos4 = field4_cos - max(field4)*cos(0:2*pi/1e2:2*pi);
figure
plot(0:2*pi/1e2:2*pi, delta_cos)
hold on
plot(0:2*pi/1e2:2*pi, delta_cos4)
clf
plot(0:2*pi/1e2:2*pi, delta_cos4-delta_cos)
plot(0:2*pi/1e2:2*pi, delta_cos4(1:51) + delta_cos4(51:101))
plot(0:2*pi/1e2:pi, delta_cos4(1:51) + delta_cos4(51:101))
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t-phi), 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
field(93)
7*(2*pi/1e2)
8*(2*pi/1e2)
phi = 0.16*pi
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t-phi), 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) maxfield*cos(t-phi), 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) maxfield*cos(t+phi), 1e-4, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) maxfield*cos(t+phi), 1e-3, options1, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
figure
plot(0:2*pi/1e2:2*pi, field)
hold on
plot(0:2*pi/1e2:2*pi, field5)
max(abs(field5-field))
save TLS2pi
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) maxfield*cos(t+phi), 1e-3, options1, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1
field_cosw = dctI(field_cos);
figure
plot(0:0.5:10, field_cosw(1:11))
plot(0:0.5:10, field_cosw(1:21))
fft_field_cos = fft(field_cos);
figure
plot(0:1:10, fft_field_cos(1:11))
fft_field_cos(1:11)
fft_field_cos = fft(field_cos(1:end-1));
plot(0:1:10, fft_field_cos(1:11))
fft_field_cos(1:11)
st_delta_cos = dstI(delta_cos(2:end-1));
figure
plot(1:1:10, st_delta_cos(1:11))
plot(1:1:10, st_delta_cos(1:10))
save TLS2pi
[P, D] = eig([0 1; 1 0])
syms phase
sigma_phi = [0 exp(1i*phase); exp(-1i*phase) 0]
eig[sigma_phi]
eig(sigma_phi)
[Pphase, Dphase] = eig(sigma_phi)
save TLS2pi
load TLS2pi
psi_max = SemiGlobalH(@(psi, t, v) H0*v - 0.5*cos(t + 0.164*pi)*sigma_x*v, @(psi1, t1, psi2, t2) 0.5*(cos(t2 + 0.164*pi) - cos(t1 + 0.164*pi)).*(sigma_x*psi1), 1, [], Edomain, [1; 0], [0  2*pi], Nt, Nt_ts, Ncheb, tol);
Edomain = [-0.2071, 1.2071]
psi_max = SemiGlobalH(@(psi, t, v) H0*v - 0.5*cos(t + 0.164*pi)*sigma_x*v, @(psi1, t1, psi2, t2) 0.5*(cos(t2 + 0.164*pi) - cos(t1 + 0.164*pi)).*(sigma_x*psi1), 1, [], Edomain, [1; 0], 0:2*pi/1e3:2*pi, Nt, Nt_ts, Ncheb, tol);
psi_max = SemiGlobalH(@(psi, t, v) H0*v - 0.5*cos(t + 0.164*pi)*sigma_x*v, @(psi1, t1, psi2, t2) 0.5*(cos(t2 + 0.164*pi) - cos(t1 + 0.164*pi)).*(sigma_x*psi1), 1, [], Edomain, [1; 0], 0:2*pi/1e3:2*pi, 1e2, 7, 7, 1e-8);
figure
plot(0:2*pi/1e3:2*pi, psi_max.*conj(pai_max)
plot(0:2*pi/1e3:2*pi, psi_max.*conj(pai_max))
plot(0:2*pi/1e3:2*pi, psi_max.*conj(psi_max))
hold on
plot(0:2*pi/1e3:2*pi, psi.*conj(psi))
plot(0:2*pi/1e2:2*pi, psi.*conj(psi))
%-- 04/04/2021 20:15 --%
load TLS2pi
whos
sigma_phi
expm(1i*sigma_phi*a)
simplify(ans)
expm(-1i*sigma_phi*a)
simplify(ans)
expm(1i*sigma_phi*a)
simplify(ans)
%-- 13/04/2021 20:17 --%
load TLS2pi
figure
plot(0:1e-2:1, fieldt)
plot(0:1e-2:1, field)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
whos
figure
plot(0:1e-2:1, field_cos)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
maxfield
ylabel('$\epsilon_{periodic}(\tau - \gamma)$', 'interpreter', 'latex')
ylabel('$\epsilon_\circ(\tau - \gamma)$', 'interpreter', 'latex')
ylabel('$\epsilon^\circ(\tau - \gamma)$', 'interpreter', 'latex')
ylabel('$\epsilon_{ext}(\tau - \gamma)$', 'interpreter', 'latex')
plot(0:1e-2:1, field_cos)
hold on
plot(0:1e-2:1, field_cos)
ylabel('$\epsilon_{ext}(\tau - \gamma)$', 'interpreter', 'latex')
hold on
plot(0:1e-2:1, cos(0:2*pi/1e2:2*pi))
plot(0:1e-2:1, cos(0:2*pi/1e2:2*pi)*maxfield)
figure
plot(0:1e-2:1, delta_cos)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
delta_cos = field_cos - maxfield*cos(0:2*pi/1e2:2*pi);
ylabel('$\epsilon_{ext}(\tau - \gamma) - \max[\epsilon(\tau)]\cos(\tau)$', 'interpreter', 'latex')
max(abs(delta_cos(1:51)+delta_cos(51:101))
max(abs(delta_cos(1:51)+delta_cos(51:101)))
max(abs(field_cos(1:51)-delta_cos(101:-1:51)))
max(abs(field_cos-delta_cos(end:-1:1)))
max(abs(field_cos-field_cos(end:-1:1)))
figure
plot(0:10, fft(field)*1/sqrt(2e2))
fftfield - fft(field)*1/sqrt(2e2);
fftfield = fft(field)*1/sqrt(2e2);
plot(0:10, fftfield(1:11))
1-J1
delta_cos(end)
delta_cos(1)
delta_cos(51)
maxfield + min(field)
fftfield = fft(field(1:end-1))*1/sqrt(2e2);
figure
plot(0:10, fftfield(1:11))
plot(0:10, fftfield(1:11).*conj(fftfield(1:11)))
max(abs(delta_cos))
max(abs(delta_cos))/maxfield
max(abs(field_cos-field_cos(end:-1:1)))
max(abs(field_cos-field_cos(end:-1:1)))/maxfield
max(abs(field_cos(1:51)+field_cos(51:101)))
fftdelta = fft(delta_cos(1:end-1))*1/sqrt(2e2);
figure
plot(0:10, fftdelta(1:11).*conj(fftdelta(1:11)))
whos
[allfield5, field5, psi5, relE5, conv5, niter5, mallniterc5, J15, maxgrad5, alpha5, invHess5] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) maxfield*cos(t+phi), 1e-3, options1, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1
figure
plot(0:1e-2:1, psi(2,:).*conj(psi(2,:)))
plot(0:1e-2:1, psi.*conj(psi))
hold on
plot(0:1e-2:1, [cos(2*pi*(0:1e-2:1)).^2; sin(2*pi*(0:1e-2:1)).^2])
plot(0:1e-2:1, psi.*conj(psi))
hold on
plot(0:1e-2:1, [cos(pi/2*(0:1e-2:1)).^2; sin(pi/2*(0:1e-2:1)).^2])
figure
plot(0:1e-2:1, psi.*conj(psi))
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
occupation = psi.*conj(psi);
0.5 - occupation(51)
0.5 - occupation(:, 51)
figure
yyaxis left
plot(phases/(2*pi)-0.5, allJmax)
plot(phases/(pi)-0.5, allJmax)
yyaxis right
plot(phases/(pi)-0.5, sin(phases))
xlabel('$\frac{\phi}{\pi}$', 'interpreter', 'latex')
yyaxis left
ylabel('$\left\vert\left<\varphi_1\vert\psi_\phi(2\pi)\right>\right\vert^2$', 'interpreter', 'latex')
ylabel('$\left\vert\left<\varphi_1\vert\psi_\phi(\tau_f)\right>\right\vert^2$', 'interpreter', 'latex')
yyaxis right
ylabel('$\cos(\phi)$', 'interpreter', 'latex')
whos
figure
plot(phases/(pi)-0.5, allJmax)
hold on
plot(phases/(pi)-0.5, allJmax_40)
plot(phases/(pi)-0.5, allJmax_45)
plot(phases/(pi)-0.5, allJmax_55)
plot(phases/(pi)-0.5, allJmax_60)
clf
plot(phases/(pi)-0.5, allJmax_40)
hold on
plot(phases/(pi)-0.5, allJmax_45)
plot(phases/(pi)-0.5, allJmax)
plot(phases/(pi)-0.5, allJmax_55)
plot(phases/(pi)-0.5, allJmax_60)
xlabel('$\frac{\phi}{\pi}$', 'interpreter', 'latex')
ylabel('$\left\vert\left<\varphi_1\vert\psi_\phi(\tau_f;\epsilon_{max})\right>\right\vert^2$', 'interpreter', 'latex')
ylabel('$\left\vert\left<\varphi_1\vert\psi(\tau_f;\phi)\right>\right\vert^2$', 'interpreter', 'latex')
ylabel('$\left\vert\left<\varphi_1\vert\psi_\phi(\tau_f;\phi,\epsilon_{max})\right>\right\vert^2$', 'interpreter', 'latex')
ylabel('$\left\vert\left<\varphi_1\vert\psi_\phi(\tau_f;\phi,\epsilon_{max}=0.5)\right>\right\vert^2$', 'interpreter', 'latex')
ylabel('$\left\vert\left<\varphi_1\vert\psi_\phi(\tau_f;\phi,0.5)\right>\right\vert^2$', 'interpreter', 'latex')
max(abs(allJmax(1:501) - allJmax(501:1001)))
max(abs(allJmax(1:500) - allJmax(501:1000)))
max(abs(allJmax(2:end) - allJmax(1000:-1:2)))
min(allJmax)
1-max(allJmax)
%-- 19/04/2021 11:40 --%
load TLS2pi
whos
expm(1i*sigma_phi*a)
simplify(ans)
U = expm(1i*sigma_phi*a)
U = simplify(U)
save TLS2pi
sigma_phi
syms Us
Us = [cos(a) 1i*exp(1i*phase)*sin(a); 1i*exp(-1i*phase)*sin(a) cos(a)]
subs(Us, a, W*tau/2)
syms W tau
subs(Us, a, W*tau/2)
Utau =  subs(Us, a, W*tau/2)
zeta = subs(sigma_phi, phase, -(2*tau + phase))
syms M1
M1 = Utau'*zeta*Utau
M1 = simplify(M1)
syms W tau phase real
M1 = simplify(M1)
M1 = Utau'*zeta*Utau
M1 = simplify(M1)
zeta
15^2
syms tau1
int(M1, tau, 0, tau)
simplify(ans)
integM1 = simplify(ans)
M1(:,1)'*M(:,1)
integM1(:,1)'*integM1(:,1)
simplify(ans)
integM12pi = subs(integM1, tau, 2*pi, W, 0.5)
integM12pi = subs(integM1, tau, 2*pi)
integM12piW = subs(integM12pi, W, 0.5)
simplify(ans)
simplify(integM12piW)
integM12piW = simplify(integM12piW)
integM12piW(:, 1)'*integM12piW(:,1)
simplify(ans)
(exp(pi/2*1i)*8i)/15 + (exp(-pi/2*3i)*8i)/15
(exp(pi/2*1i)*8i)/15 + (exp(-3*pi/2*3i)*8i)/15
exp(-pi/2*3i)
-136*15/(8*225)
16/225
8*15
U2piW =  subs(Utau, tau, 2*pi)
U2piW =  subs(U2piW, W, 0.5)
U1= 1i/4*U2piW*integM12piW
U1=simplify(U1)
U1(:,1)'*U1(:,1)
simplify(ans)
-30*17/450
17/450
1/30
17/450-1/30
save TLS2pi
136/225-8/15
ans/4
7.1111e-02/16
U1gen = 1i/4*Utau*integM1
U1gen = simplify(U1gen)
U1W = subs(U1gen, W, 0.5)
U1W = simplify(U1W)
M2 = Utau'*zeta*U1W
M2 = Utau'*zeta*U1gen
M2 = simplify(M2)
M2W = subs(M2W, W, 0.5)
M2W = subs(M2, W, 0.5)
M2W = simplify(M2W)
U1gen = 1i*W/2*Utau*integM1
U1gen = simplify(U1gen)
U1W = subs(U1gen, W, 0.5)
U1W = simplify(U1W)
M2 = Utau'*zeta*U1gen
M2 = simplify(M2)
integM2 = int(M2, tau, 0, tau)
integM2 = simplify(integM2)
U2gen= 1i*W/2*Utau*integM2
U2gen=simplify(U2gen)
U2W = subs(U2gen, W, 0.5)
U2W = simplify(U2W)
U2W2pi = subs(U2W, tau, 2*pi)
U2W2pi = simplify(U2W2pi)
U2 = simplify(U2W2pi)
U12 = U1+U2
U12(:,1)'*U12(:,1)
U12sqn = simplify(U12(:,1)'*U12(:,1))
DU12sqn = diff(U12sqn)
D2U12sqn = diff(DU12sqn)
solve(DU12sqn == 0,phase)
vpasolve(DU12sqn == 0,phase)
-acos(root(z^3 + 30*z^2 + (3047*z)/5 - 330, z, 1))/2
syms(z)
syms z
-acos(root(z^3 + 30*z^2 + (3047*z)/5 - 330, z, 1))/2
root(z^3 + 30*z^2 + (3047*z)/5 - 330, z, 1)
subs(D2U12sqn, 0)
subs(D2U12sqn, pi/2)
subs(D2U12sqn, -pi/2)
subs(DU12sqn, -pi/2)
subs(DU12sqn, pi/2)
subs(DU12sqn, 0)
solve(DU12sqn == 0,phase, 0.16*pi)
solve(DU12sqn == 0,phase, [0.16*pi, 0.17*pi)
solve(DU12sqn == 0,phase, [0.16*pi, 0.17*pi])
solve(DU12sqn == 0,phase, 0.164*pi)
solve(z^3 + 30*z^2 + (3047*z)/5 - 330==0, z)
vpasolve(z^3 + 30*z^2 + (3047*z)/5 - 330==0, z)
-acos(0.52757327196539757966794748921944)/2
subs(DU12sqn, ans)
subs(ans)
simplify(ans)
vpa(ans)
-5.0753e-01/pi
subs(DU12sqn, ans)
subs(D2U12sqn, -1.6155e-01)
vpa(ans)
subs(D2U12sqn, -1.6155e-01*pi)
vpa(ans)
extremal_phase = solve(DU12sqn == 0,phase)
vpa(extremal_phase)
extremal_phase_num = vpa(extremal_phase)
subs(DU12sqn, extremal_phase_num)
subs(D2U12sqn, extremal_phase_num)
subs(U12sqn, extremal_phase_num)
vpn(ans)
vpa(ans)
solve(DU12sqn == 0,phase, 'MaxDegree', 3)
vpa(ans)
subs(DU12sqn, -pi/2)
subs(U12sqn, [-pi/2; pi/2])
vpa(ans)
subs(DU12sqn, [-pi/2; pi/2])
subs(D2U12sqn, [-pi/2; pi/2])
subs(D2U12sqn, extremal_phase_num)
doc subs
sampU12sqn = vps(subs(U12sqn, [-pi/2:pi/1e3:3*pi/2]));
sampU12sqn = vpa(subs(U12sqn, [-pi/2:pi/1e3:3*pi/2]));
figure
plot(-pi/2:pi/1e3:3*pi/2, 1-sampU12sqn)
plot(-1/2:1/1e3:3/2, 1-sampU12sqn)
hold on
plot(phases/(pi)-0.5, allJmax)
save TLS2pi
U12oc0 = simplify(U12(1,1)*conj(U12(1,1)))
DU12oc0 = diff(U12oc0)
D2U12oc0 = diff(DU12oc0)
extremal_phase1 = solve(DU12oc0 == 0,phase)
extremal_phase1_num = vpa(extremal_phase1)
extremal_phase1_num/pi
vpa(extremal_phase1_num/pi)
subs(D2U12oc0, extremal_phase1_num)
subs(DU12oc0, extremal_phase1_num)
subs(DU12oc0, [-pi/2])
subs(DU12oc0, [pi/2])
sampU12oc0 = vpa(subs(U12oc0, [-pi/2:pi/1e3:3*pi/2]));
plot(-1/2:1/1e3:3/2, 1-sampU12oc0)
save TLS2pi
%-- 22/04/2021 14:46 --%
U1
load TLS2pi
U1
U1gen
U1gen(1,1)
U2gen(1,1)
U2
U2W
subs(U2W, tau, 0)
subs(U2gen, tau, 0)
subs(U1gen, tau, 0)
simplify(ans)
subs(U2gen, W, pi/tau)
simplify(ans)
subs(ans, Inf)
max(abs(allJmax))
max(abs(allJmax6))
whos
max(abs(allJmax_60))
(1-cos(0.6*2*pi))/2
max(abs(allJmax_40))
(1-cos(0.4*2*pi))/2
max(abs(allJmax_45))
(1-cos(0.45*2*pi))/2
max(abs(allJmax_55))
(1-cos(0.55*2*pi))/2
M2
zeta
integM2
doc int
U1gen
integM1
int(M2, tau, 0, tauf)
M2tau1 = subs(M2, tau1)
int(M2tau1, tau1, 0, tau)
integM2b = simplify(ans)
integM2b-integM2
clear integM2b
fftfield(1:10)*2*pi/sqrt(2e2)
whos
fftfield = fft(field(1:end-1))*1/sqrt(2e2);
ftfield = fft(field(1:end-1))*2*pi/2e2;
ftfield(1:10)
ftfield(1:10).*conj(ftfield(1:10))
fftfield(1:10).*conj(fftfield(1:10))
abs(ftfield(1:10))
doc fft
ftfield(1001)
length(field)
ftfield(101)
ftfield(51)
ftfield(52)
ftfield(53)
ftfield(50)
figure
plot(0:100, ftfield.*conj(ftfield))
size(ftfield)
plot(0:99, ftfield.*conj(ftfield))
ftfield(100)
ftfield(2)
abs(ans)
log(ftfield(2)/abs(ftfield(2))
log(ftfield(2)/abs(ftfield(2)))
ans/pi
2*real(ftfield(2))/(2*pi)
2*real(ftfield(2))
2*abs(ftfield(2))/(2*pi)
2*abs(ftfield(2))/(pi)
ftfield = fft(field(1:end-1))*sqrt(2*pi)/2e2;
size(field)
ftfield = fft(field(1:end-1))*sqrt(2*pi)/100;
2*real(ftfield(2))/sqrt(2*pi)
2*abs(ftfield(2))/sqrt(2*pi)
2*abs(fftfield(2))/100
fftfield = fft(field(1:end-1));
2*abs(fftfield(2))/100
Afield_w1 = 2*abs(fftfield(2))/100
save TLS2pi
allJmax_45 = JmaxTLSphases(0.45, 1e3, [-0.2071, 1.2071], 100, 7, 7, 1e-8);
allJmax_Aw1 = JmaxTLSphases(Afield_w1, 1e3, [-0.2071, 1.2071], 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax_Aw1)
1-max(allJmax_Aw1)
allJmax_Aw1fine = JmaxTLSphasesD(Afield_w1, [331, 333]*2*pi/1e3, 1001, [-0.2071, 1.2071], 100, 7, 7, 1e-8);
figure
plot(phases_fine, allJmax_Aw1fine)
1-max(allJmax_Aw1fine)
figure
plot(0:10, log(abs(ftfield(1:11))))
plot(0:15, log(abs(ftfield(1:16))))
plot(0:10, log(abs(ftfield(1:11))))
log(fftfield(4)/abs(fftfield(4)))
log(fftfield(4)/abs(fftfield(4)))/1i
log(fftfield(2)/abs(fftfield(2)))/1i
log(fftfield(4)/abs(fftfield(4)))/(1i*pi)
log(fftfield(6)/abs(fftfield(6)))/(1i*pi)
fftfield_minus1 = fftfield;
fftfield_minus1([2 100]) = 0;
field_minus1 = ifft(fftfield_minus1);
figure
plot(0:100,fftfield)
plot(0:99,fftfield)
plot(0:99,abs(fftfield))
plot(0:99,abs(fftfield_minus1))
plot(0:99,abs(fftfield))
plot(0:99,abs(fftfield_minus1))
plot(0:2*pi/1e2:2*pi, field_minus1)
plot(0:2*pi/1e2:99/100*2*pi, field_minus1)
figure
plot(0:2*pi/1e2:99/100*2*pi, field_minus1)
phase_w1 = log(fftfield(6)/abs(fftfield(6)))/(-1i)
phase_w1 = log(fftfield(6)/abs(fftfield(6)))/(1i)
phase_w1 = real(log(fftfield(2)/abs(fftfield(2)))/(1i))
phase_w3 = real(log(fftfield(4)/abs(fftfield(4)))/(1i))
phase_w5 = real(log(fftfield(6)/abs(fftfield(6)))/(1i))
phase_w3/pi
phase_w5/pi
abs(fftfield(1:10))
abs(fftfield(4))/abs(fftfield(2))
plot(phases_fine/pi - 0.5, allJmax_Aw1fine)
[maxJmaxAw1fine, imaxJmaxAw1fine] = max(allJmax_Aw1fine)
1-maxJmaxAw1fine
allJmax_Aw1fine(615)
1-allJmax_Aw1fine(615)
1-allJmax_Aw1fine(614)
1-allJmax_Aw1fine(616)
allJmax_Aw1(26)
allJmax_Aw1(251)
allJmax_Aw1(331)
allJmax_Aw1(333)
-0.5+330*1e-3+
-0.5+330*1e-3
-0.5+330*1e-3+614*2e-6
-0.5+331*1e-3+614*2e-6
-0.5+331*2e-3+614*2e-6
-0.5+331*2e-3+614*4e-6
xlabel('$\nu$', 'interpreter', 'latex')
ylabel('$\left\vert\bar\epsilon(\nu)\right\vert$', 'interpreter', 'latex')
plot(0:10, log10(abs(ftfield(1:11))))
xlabel('$\nu$', 'interpreter', 'latex')
ylabel('$\left\vert\bar\epsilon(\nu)\right\vert$', 'interpreter', 'latex')
phase_w1/pi
1-allJmax_Aw1fine(19)
0.00018/4e-6
1-allJmax_Aw1fine(46)
ifftfield = ifft(field(1:end-1));
abs(ifftfield(1:10))
2*abs(ifftfield(1:10))
[ifftfield(1:10); fftfield(1:10)*100]
[ifftfield(1:10); fftfield(1:10)/100]
plot(0:10, log10(2*abs(ifftfield(1:11))))
real(log(ifftfield(2)/abs(ifftfield(2)))/(-1i))
xlabel('$\nu$', 'interpreter', 'latex')
ylabel('$\log_{10}\left(2\left\vert\bar\epsilon(\nu)\right\vert\right)$', 'interpreter', 'latex')
figure
plot(0:1e-2:1, field)
hold on
plot(0:1e-2:1, 2*ifftfield(2)*cos((0:pi/100:pi)) + phase_w1)
plot(0:1e-2:1, 2*abs(ifftfield(2))*cos((0:2*pi/100:2*pi)) + phase_w1)
abs(ifftfield(2))
plot(0:1e-2:1, 2*abs(ifftfield(2))*cos((0:2*pi/100:2*pi) + phase_w1))
whos
figure
plot(phases/(pi)-0.5, allJmax)
hold on
plot(-1/2:1/1e3:3/2, 1-sampU12oc0)
1-max(sampU12oc0)
1-min(sampU12oc0)
save TLS2pi
%-- 03/05/2021 11:03 --%
load TLS2pi
U1
U1gen
U2
U1
U2gen
integM1
U
U1gen = 1i/4*
Utau*integM1
zeta
logm(zeta)
simplify(ans)
expm(pi*1i/2*eye(2) - zeta)
simplify(ans)
expm(pi*1i/2*(eye(2) - zeta))
ans-zeta
whos
U1
1/15
U2
pi/120
(1/15)^2
(pi/120 + 1/315 + 2*i/15)
abs(pi/120 + 1/315 + 2*i/15)^2
abs(pi/120 + 2/315 + 2*i/15)^2
abs(pi/120 + 2/315 + 2/15)^2
U12oc0
pi^2/14400 + 1936/99225
max(sampU12oc0)
6704/99225
352/4725
2/15
(2/15)^2
3/15
U1phase0 = subs(U1gen, phase, 0)
U1phase0 = simplify(subs(U1gen, phase, 0))
U1Wphase0 = simplify(subs(U1phase0, W, 0.5))
U2phase0 = simplify(subs(U2gen, phase, 0))
U2Wphase0 = simplify(subs(U2phase0, W, 0.5))
sqnU1Wphase0 = simplify(U1Wphase0(:,1)'*U1Wphase0(:,1))
samp_sqnU1Wphase0 = vpa(subs(sqnU1Wphase0, 0:pi/1e3:2*pi));
sqnU2Wphase0 = simplify(U2Wphase0(:,1)'*U2Wphase0(:,1))
samp_sqnU2Wphase0 = vpa(subs(sqnU2Wphase0, 0:pi/1e3:2*pi));
samp_sqnU2Wphase0 = vpa(simplify(subs(sqnU2Wphase0, pi/1e3:pi/1e3:2*pi)));
samp_sqnU2Wphase0 = vpa((subs(sqnU2Wphase0, pi/1e3:pi/1e3:2*pi)));
figure
plot(0:pi/1e3:2*pi, samp_sqnU1Wphase0)
plot(0:pi/1e3:2*pi, imag(samp_sqnU1Wphase0))
plot(0:pi/1e3:2*pi, real(samp_sqnU1Wphase0))
hold on
plot(pi/1e3:pi/1e3:2*pi, real(samp_sqnU2Wphase0))
whos
sigma_x*zeta*sigma_x
zeta
sigma_x*conj(zeta)*sigma_x
sigma_x*conj(zeta)*sigma_x-zeta
sigma_x*zeta
zeta*sigma_x
zeta
sigma_x*[1 2; 3 4]
[1 2; 3 4]*sigma_x
U
save TLS2pi
[allfield4pi, field4pi, psi4pi, relE4pi, conv4pi, niter4pi, mallniterc4pi, J14pi, maxgrad4pi, alpha4pi, invHess4pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*cos(t), 1e-3, options, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfield4pi, field4pi, psi4pi, relE4pi, conv4pi, niter4pi, mallniterc4pi, J14pi, maxgrad4pi, alpha4pi, invHess4pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*cos(t), 5e-4, options, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
J14pi
1-J14pi
figure
plot(0:pi/50:4*pi, field4pi)
options4pi = optionsOCqn(1e-4,1e4);
options4pi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 0.4)
[allfield4pi, field4pi, psi4pi, relE4pi, conv4pi, niter4pi, mallniterc4pi, J14pi, maxgrad4pi, alpha4pi, invHess4pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*cos(t), 5e-4, options4pi, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
options4pi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 0.5)
[allfield4pi, field4pi, psi4pi, relE4pi, conv4pi, niter4pi, mallniterc4pi, J14pi, maxgrad4pi, alpha4pi, invHess4pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*cos(t), 5e-4, options4pi, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J14pi
figure
plot(0:pi/50:4*pi, field4pi)
ifftfield4pi = ifft(field4pi(1:end-1));
figure
plot(0:10, log10(2*abs(ifftfield4pi(1:11))))
plot(0:0.5:10, log10(2*abs(ifftfield4pi(1:21))))
[allfield4pi, field4pi, psi4pi, relE4pi, conv4pi, niter4pi, mallniterc4pi, J14pi, maxgrad4pi, alpha4pi, invHess4pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*cos(t), 1e-4, options4pi, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
plot(0:pi/50:4*pi, field4pi)
1-J14pi
ifftfield4pi = ifft(field4pi(1:end-1));
plot(0:0.5:10, log10(2*abs(ifftfield4pi(1:21))))
dctIfield4pi = dctI(field4pi);
figure
plot(0:0.25:10, field4pi(1:41))
[allfield4pi1, field4pi1, psi4pi1, relE4pi1, conv4pi1, niter4pi1, mallniterc4pi1, J14pi1, maxgrad4pi1, alpha4pi1, invHess4pi1] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*sin(t), 1e-4, options4pi, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
options4pi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 1)
[allfield4pi1, field4pi1, psi4pi1, relE4pi1, conv4pi1, niter4pi1, mallniterc4pi1, J14pi1, maxgrad4pi1, alpha4pi1, invHess4pi1] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 0.25*sin(t), 1e-4, options4pi, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J14pi1
1-J14pi
figure
conv4pi(end)
conv4pi1(end)
figure
plot(0:0.25:10, ifftfield4pi(1:41))
plot(0:0.25:10, log10(2*abs(ifftfield4pi(1:41))))
plot(0:0.25:10, log10(abs(dctIfield4pi(1:41))))
plot(0:pi/50:4*pi, field4pi1)
conv4pi1(end)-conv4pi(end)
ifftfield4pi(1:11)
2*ifftfield4pi(1:11)
allJmax4pi = JmaxTLSphasesT(0.25, 0:pi/500:(2*pi*999/1e3), [-6 6], 4*pi, 100, 7, 7, 1e-8);
size(allJmax4pi)
size(allJmax)
figure
plot(phases/(pi)-0.5, allJmax4pi)
allJmax3pi = JmaxTLSphasesT(1/3, 0:pi/500:(2*pi*999/1e3), [-6 6], 3*pi, 75, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax3pi)
allJmax2p5pi = JmaxTLSphasesT(1/2.5, 0:pi/500:(2*pi*999/1e3), [-6 6], 2.5*pi, 75, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax2p5pi)
allJmaxpi = JmaxTLSphasesT(1, 0:pi/500:(2*pi*999/1e3), [-6 6], pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmaxpi)
save TLS2pi
allJmax3p5pi = JmaxTLSphasesT(1/3.5, 0:pi/500:(2*pi*999/1e3), [-6 6], 3.5*pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax3p5pi)
allJmax05pi = JmaxTLSphasesT(2, 0:pi/500:(2*pi*999/1e3), [-6 6], pi/2, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax05pi)
%-- 11/05/2021 15:36 --%
allJmax3p5pi = JmaxTLSphasesT(1/3.5, 0:pi/500:(2*pi*999/1e3), [-6 6], 3.5*pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax3p5pi)
allJmax05pi = JmaxTLSphasesT(2, 0:pi/500:(2*pi*999/1e3), [-6 6], pi/2, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax05pi)
load TLS2pi
allJmax3p5pi = JmaxTLSphasesT(1/3.5, 0:pi/500:(2*pi*999/1e3), [-6 6], 3.5*pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax3p5pi)
allJmax05pi = JmaxTLSphasesT(2, 0:pi/500:(2*pi*999/1e3), [-6 6], pi/2, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax05pi)
allJmax8pi = JmaxTLSphasesT(1/8, 0:pi/500:(2*pi*999/1e3), [-6 6], 8*pi, 200, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax8pi)
allJmax8pi_fine = JmaxTLSphasesT(1/8, 0.165*pi:pi/1e5:0.168, [-6 6], 8*pi, 200, 7, 7, 1e-8);
size(allJmax8pi_fine)
allJmax8pi_fine = JmaxTLSphasesT(1/8, 0.165*pi:pi/1e5:0.168*pi, [-6 6], 8*pi, 200, 7, 7, 1e-8);
size(allJmax8pi_fine)
figure
plot(0.165:1e-55:0.168, allJmax8pi_fine)
plot(0.165:1e-5:0.168, allJmax8pi_fine)
allJmax8pi_fine = JmaxTLSphasesT(1/8, 0.165*pi:pi/1e5:0.17*pi, [-6 6], 8*pi, 200, 7, 7, 1e-8);
plot(0.165:1e-5:0.17, allJmax8pi_fine)
allJmax8pi_fine = JmaxTLSphasesT(1/8, pi/2+(0.165*pi:pi/1e5:0.168*pi), [-6 6], 8*pi, 200, 7, 7, 1e-8);
plot(0.165:1e-5:0.168, allJmax8pi_fine)
allJmax16pi_fine = JmaxTLSphasesT(1/16, 0.165*pi:pi/1e5:0.17*pi, [-6 6], 16*pi, 400, 7, 7, 1e-8);
figure
plot(0.165:1e-5:0.168, allJmax16pi_fine)
plot(0.165:1e-5:0.17, allJmax16pi_fine)
allJmax16pi_fine = JmaxTLSphasesT(1/16, pi/2+(0.165*pi:pi/1e5:0.168*pi), [-6 6], 16*pi, 200, 7, 7, 1e-8);
allJmax16pi_fine = JmaxTLSphasesT(1/16, pi/2+(0.165*pi:pi/1e5:0.168*pi), [-6 6], 16*pi, 100, 7, 7, 1e-8);
plot(0.165:1e-5:0.17, allJmax16pi_fine)
plot(0.165:1e-5:0.168, allJmax16pi_fine)
[max16fine, imax16fine] = max(allJmax16pi_fine)
0.165*163*1e-5
0.165+163*1e-5
allJmax32pi_fine = JmaxTLSphasesT(1/32, pi/2+(0.165*pi:pi/1e5:0.168*pi), [-6 6], 32*pi, 100, 7, 7, 1e-8);
allJmax32pi_fine = JmaxTLSphasesT(1/32, pi/2+(0.165*pi:pi/1e5:0.168*pi), [-6 6], 32*pi, 100, 7, 9, 1e-8);
allJmax32pi_fine = JmaxTLSphasesT(1/32, pi/2+(0.165*pi:pi/1e5:0.168*pi), [-6 6], 32*pi, 200, 7, 7, 1e-8);
figure
plot(0.165:1e-5:0.168, allJmax32pi_fine)
[max32fine, imax32fine] = max(allJmax32pi_fine)
0.165+166*1e-5
U1gen
U1genfWopt = subs(U1gen, W, pi/tauf)
syms tauf
U1genfWopt = subs(U1gen, tau, tauf)
U1genfWopt = subs(U1genfWopt, W, pi/tauf)
U1genfWopt = simplify(U1genfWopt)
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, options4pi, pi, 2*pi/25, 7, 7, 1e-4, 1e4);
optionspi = optionsOCqn(1e-4,1e4);
optionspi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 4)
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, optionspi, pi, 2*pi/25, 7, 7, 1e-4, 1e4);
1-J1pi
figure
plot(0:pi/50:pi, fieldpi)
size(fieldpi)
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, optionspi, pi, pi/50, 7, 7, 1e-4, 1e4);
optionspi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 8)
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, optionspi, pi, pi/50, 7, 7, 1e-4, 1e4);
optionspi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 16)
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, optionspi, pi, pi/50, 7, 7, 1e-4, 1e4);
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, options, pi, pi/50, 7, 7, 1e-4, 1e4);
optionspi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 32)
[allfieldpi, fieldpi, psipi, relEpi, convpi, niterpi, mallnitercpi, J1pi, maxgradpi, alphapi, invHesspi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t), 1e-4, optionspi, pi, pi/50, 7, 7, 1e-4, 1e4);
1-J1pi
figure
plot(0:pi/50:pi, fieldpi)
figure
plot(phases/(pi)-0.5, allJmaxpi)
[allfieldpi1, fieldpi1, psipi1, relEpi1, convpi1, niterpi1, mallnitercpi1, J1pi1, maxgradpi1, alphapi1, invHesspi1] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) cos(t + 0.158*pi), 1e-4, optionspi, pi, pi/50, 7, 7, 1e-4, 1e4);
1-J1pi1
1-J1pi
convpi(end)
convpi1(end)
convpi1(end)-convpi(end)
figure
plot(0:pi/50:pi, fieldpi1)
fieldpi1ext = [fieldpi1(1:end-1), -fieldpi1(1:end-1)];
figure
plot(0:pi/50:1.98*pi, fieldpi1ext)
ifftfieldpi1ext = ifft(fieldpi1ext);
figure
plot(0:10, log10(2*abs(ifftfieldpi1ext(1:11))))
ifftfieldpi1ext(1:10)
plot(0:10, log10(2*abs(ifftfieldpi1ext(1:11))))
abs(ifftfieldpi1ext(1:11))
2*ans
plot(1:2:21, log10(2*abs(ifftfieldpi1ext(2:2:22))))
[allfield02pi, field02pi, psi02pi, relE02pi, conv02pi, niter02pi, mallniterc02pi, J102pi, maxgrad02pi, alpha02pi, invHess02pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 5*cos(t), 1e-3, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
1-J102pi
[allfield02pi, field02pi, psi02pi, relE02pi, conv02pi, niter02pi, mallniterc02pi, J102pi, maxgrad02pi, alpha02pi, invHess02pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 5*cos(t), 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pi1, field02pi1, psi02pi1, relE02pi1, conv02pi1, niter02pi1, mallniterc02pi1, J102pi1, maxgrad02pi1, alpha02pi1, invHess02pi1] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield02pi, 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pi, field02pi, psi02pi, relE02pi, conv02pi, niter02pi, mallniterc02pi, J102pi, maxgrad02pi, alpha02pi, invHess02pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 5*cos(t), 1e-3, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pi1, field02pi1, psi02pi1, relE02pi1, conv02pi1, niter02pi1, mallniterc02pi1, J102pi1, maxgrad02pi1, alpha02pi1, invHess02pi1] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield02pi, 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pi1, field02pi1, psi02pi1, relE02pi1, conv02pi1, niter02pi1, mallniterc02pi1, J102pi1, maxgrad02pi1, alpha02pi1, invHess02pi1] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield02pi.', 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
1-J102pi1
figure
plot(0:pi/200:pi/5, field02pi)
hold on
plot(0:pi/200:pi/5, field02pi1)
[allfield02pi2, field02pi2, psi02pi2, relE02pi2, conv02pi2, niter02pi2, mallniterc02pi2, J102pi2, maxgrad02pi2, alpha02pi2, invHess02pi2] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield02pi1.', 1e-5, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
1-J102pi2
plot(0:pi/200:pi/5, field02pi2)
figure
plot(0:pi/200:pi/5, field02pi2)
[allfield02pi2, field02pi2, psi02pi2, relE02pi2, conv02pi2, niter02pi2, mallniterc02pi2, J102pi2, maxgrad02pi2, alpha02pi2, invHess02pi2] = OCqn([1;0], [0;1], H0, [-6 6], mu, allfield02pi1.', 1e-6, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pi2, field02pi2, psi02pi2, relE02pi2, conv02pi2, niter02pi2, mallniterc02pi2, J102pi2, maxgrad02pi2, alpha02pi2, invHess02pi2] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 5*cos(5*(t-0.141)), 1e-5, optionspi, 0.2*pi, pi/400, 7, 7, 1e-4, 1e4);
1-J102pi2
figure
plot(0:pi/200:pi/5, field02pi2)
plot(0:pi/400:pi/5, field02pi2)
field02pi2ext = [field02pi2(1:end-1), -field02pi2(1:end-1)];
figure
plot(0:pi/400:(0.4-1/400)*pi, field02pi2ext)
ifftfieldpi1ext = ifft(field02pi2ext);
figure
ifftfieldpi1ext(1:10)
2*abs(ifftfieldpi1ext(1:10))
arg(ifftfieldpi1ext(1:10))
angle(ifftfieldpi1ext(1:10))
angle(ifftfieldpi1ext(1:10))/pi
angle(ifftfieldpi1ext(1:10))*0.2/pi
angle(ifftfieldpi1ext(1:10))*0.4/pi
log(ifftfieldpi1ext(1:10))/(-1i)
log(ifftfieldpi1ext(1:10)./abs(ifftfieldpi1ext(1:10)))/(-1i)
log(ifftfieldpi1ext(1:10)./abs(ifftfieldpi1ext(1:10)))/(-1i)/pi
save TLS2pi
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-5, optionspi, 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-3, optionspi, 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
options01pi = optionsOCqn(1e-4,1e4);
options01pi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 50)
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-3, optionspi, 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(t), 1e-3, optionspi, 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-3, options01pi, 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
1-J101pi
figure
plot(0:pi/400:pi/10, field01pi)
[allfield015pi, field015pi, psi015pi, relE015pi, conv015pi, niter015pi, mallniterc015pi, J1015pi, maxgrad015pi, alpha015pi, invHess015pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 1/0.15*cos(t/0.15), 1e-3, options01pi, 0.15*pi, 0.0015*pi, 7, 7, 1e-4, 1e4);
1-J1015pi
[allfield015pia, field015pia, psi015pia, relE015pia, conv015pia, niter015pia, mallniterc015pia, J1015pia, maxgrad015pia, alpha015pia, invHess015pia] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield015pi, 1e-3, options01pi, 0.15*pi, 0.0015*pi, 7, 7, 1e-4, 1e4);
[allfield015pia, field015pia, psi015pia, relE015pia, conv015pia, niter015pia, mallniterc015pia, J1015pia, maxgrad015pia, alpha015pia, invHess015pia] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield015pi.', 1e-3, options01pi, 0.15*pi, 0.0015*pi, 7, 7, 1e-4, 1e4);
figure
plot(0:0.0015*pi:0.15*pi, field015pi)
hold on
plot(0:0.0015*pi:0.15*pi, field015pia)
1/0.15
[allfield01pia, field01pia, psi01pia, relE01pia, conv01pia, niter01pia, mallniterc01pia, J101pia, maxgrad01pia, alpha01pia, invHess01pia] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield015pi, 1e-3, options01pi, 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
[allfield01pia, field01pia, psi01pia, relE01pia, conv01pia, niter01pia, mallniterc01pia, J101pia, maxgrad01pia, alpha01pia, invHess01pia] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield015pi, 1e-3, options01pi.', 0.1*pi, pi/400, 7, 7, 1e-4, 1e4);
[allfield015pib, field015pib, psi015pib, relE015pib, conv015pib, niter015pib, mallniterc015pib, J1015pib, maxgrad015pib, alpha015pib, invHess015pib] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield015pia.', 1e-4, options01pi, 0.15*pi, 0.0015*pi, 7, 7, 1e-4, 1e4);
1-J1015pia
1-J1015pib
1-J1015pi
figure
plot(0:0.0015*pi:0.15*pi, field015pib)
[allfield015pic, field015pic, psi015pic, relE015pic, conv015pic, niter015pic, mallniterc015pic, J1015pic, maxgrad015pic, alpha015pic, invHess015pic] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield015pib.', 1e-5, options01pi, 0.15*pi, 0.0015*pi, 7, 7, 1e-4, 1e4);
1-J1015pic
figure
plot(0:0.0015*pi:0.15*pi, field015pic)
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-2, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-2, options01pi, 0.1*pi, 5e-4*pi, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-1, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
options01pi.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 100)
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-1, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-4, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pi, field01pi, psi01pi, relE01pi, conv01pi, niter01pi, mallniterc01pi, J101pi, maxgrad01pi, alpha01pi, invHess01pi] = OCqn([1;0], [0;1], H0, [-15 15], mu, @(t) 10*cos(10*t), 1e-3, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
1-J101pi
figure
plot(0:1e-3*pi:0.1*pi, field01pi)
[allfield01pi1, field01pi1, psi01pi1, relE01pi1, conv01pi1, niter01pi1, mallniterc01pi1, J101pi1, maxgrad01pi1, alpha01pi1, invHess01pi1] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield01pi, 1e-5, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pi1, field01pi1, psi01pi1, relE01pi1, conv01pi1, niter01pi1, mallniterc01pi1, J101pi1, maxgrad01pi1, alpha01pi1, invHess01pi1] = OCqn([1;0], [0;1], H0, [-15 15], mu, allfield01pi.', 1e-5, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
1-J101pi1
figure
plot(0:1e-3*pi:0.1*pi, field01pi1)
save TLS2pi
field01pi1ext = [field02pi1(1:end-1), -field01pi1(1:end-1)];
ifftfield01pi1ext = ifft(field01pi1ext);
ifftfieldpi1ext = ifft(fieldpi1ext);
ifftfield02pi2ext = ifft(field02pi2ext);
log(ifftfield02pi2ext(1:10)./abs(ifftfield02pi2ext(1:10)))/(-1i)/pi
log(ifftfieldpi1ext(1:10)./abs(ifftfieldpi1ext(1:10)))/(-1i)/pi
log(ifftfield02pi2ext(1:10)./abs(ifftfield02pi2ext(1:10)))*0.4/(-1i)/pi
figure
2*abs(ifftfieldpi1ext(1:10))
2*abs(ifftfield02pi2ext(1:10))
2*abs(ifftfield01pi1ext(1:10))
field01pi1ext = [field01pi1(1:end-1), -field01pi1(1:end-1)];
ifftfield02pi2ext = ifft(field02pi2ext);
log(ifftfield01pi1ext(1:10)./abs(ifftfield01pi1ext(1:10)))*0.05/(-1i)/pi
log(ifftfield01pi1ext(1:10)./abs(ifftfield01pi1ext(1:10)))/(-1i)/pi
log(ifftfield01pi1ext(1:10)./abs(ifftfield01pi1ext(1:10)))*20/(-1i)/pi
field01pi1ext = [field01pi1(1:end-1), -field01pi1(1:end-1)];
plot(0:1e-3*pi:0.099*pi, field01pi1ext)
size(field01pi1ext)
plot(0:1e-3*pi:0.199*pi, field01pi1ext)
ifftfield01pi1ext = ifft(field01pi1ext);
ifftfield01pi1ext(1:11)
2*abs(ifftfield01pi1ext(1:11))
log(ifftfield01pi1ext(1:10)./abs(ifftfield01pi1ext(1:10)))/(-1i)/pi
log(ifftfieldpi1ext(1:10)./abs(ifftfieldpi1ext(1:10)))/(-1i)/pi
log(ifftfield02pi2ext(1:10)./abs(ifftfield02pi2ext(1:10)))/(-1i)/pi
0.1 + log(ifftfield02pi2ext(1:10)./abs(ifftfield02pi2ext(1:10)))/(-1i)/pi
[allfield2p5pi, field2p5pi, psi2p5pi, relE2p5pi, conv2p5pi, niter2p5pi, mallniterc2p5pi, J12p5pi, maxgrad2p5pi, alpha2p5pi, invHess2p5pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 1/2.5*cos(t), 1e-3, options, 2.5*pi, pi/50, 7, 7, 1e-4, 1e4);
[allfield2p5pi, field2p5pi, psi2p5pi, relE2p5pi, conv2p5pi, niter2p5pi, mallniterc2p5pi, J12p5pi, maxgrad2p5pi, alpha2p5pi, invHess2p5pi] = OCqn([1;0], [0;1], H0, [-6 6], mu, @(t) 1/2.5*cos(t), 1e-3, options, 2.5*pi, pi/100, 7, 7, 1e-4, 1e4);
1-J12p5pi
figure
plot(0:pi/100:2.5*pi, field2p5pi)
0.9739/pi
5.9/(2*pi)
allJmax2p01pi = JmaxTLSphasesT(1/2.01, 0:pi/500:(2*pi*999/1e3), [-6 6], 2.01*pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax2p01pi)
allJmax2p1pi = JmaxTLSphasesT(1/2.1, 0:pi/500:(2*pi*999/1e3), [-6 6], 2.1*pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmax2p1pi)
save TLS2pi
muxz = [1 1; 1 -1]
[allfieldxz, fieldxz, psixz, relExz, convxz, niterxz, mallnitercxz, J1xz, maxgradxz, alphaxz, invHessxz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 0.5*cos(t), 1e-3, options, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
[allfieldxz, fieldxz, psixz, relExz, convxz, niterxz, mallnitercxz, J1xz, maxgradxz, alphaxz, invHessxz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 0.5*cos(t), 1e-3, options, 2*pi, pi/1e2, 7, 7, 1e-4, 1e4);
[allfieldxz, fieldxz, psixz, relExz, convxz, niterxz, mallnitercxz, J1xz, maxgradxz, alphaxz, invHessxz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 0.5*cos(t), 1e-3, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1xz
figure
plot(0:pi/100:2*pi, fieldxz)
plot(0:pi/50:2*pi, fieldxz)
[allfieldxz, fieldxz, psixz, relExz, convxz, niterxz, mallnitercxz, J1xz, maxgradxz, alphaxz, invHessxz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 0.5*cos(t), 1e-4, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1xz
plot(0:pi/50:2*pi, fieldxz)
mean(fieldxz)
plot(0:pi/50:2*pi, fieldxz-mean(fieldxz))
plot(0:pi/50:2*pi, fieldxz)
[allfieldxz4, fieldxz4, psixz4, relExz4, convxz4, niterxz4, mallnitercxz4, J1xz4, maxgradxz4, alphaxz4, invHessxz4] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 0.25*cos(t), 1e-4, optionspi, 4*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1xz4
figure
plot(0:pi/50:4*pi, fieldxz4)
allJmaxxz = JmaxTLSphases_mu(muxz, 1/2, 0:pi/500:(2*pi*999/1e3), [-6 6], 2*pi, 100, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmaxxz)
allJmaxxz4 = JmaxTLSphases_mu(muxz, 0.25, 0:pi/500:(2*pi*999/1e3), [-6 6], 4*pi, 200, 7, 7, 1e-8);
figure
plot(phases/(pi)-0.5, allJmaxxz4)
save TLS2pi
U1genfWopt
U1
save TLS2pi
[allfield01pixz, field01pixz, psi01pixz, relE01pixz, conv01pixz, niter01pixz, mallniterc01pixz, J101pixz, maxgrad01pixz, alpha01pixz, invHess01pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 10*cos(10*t), 1e-3, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
1-J101pixz
[allfield01pixza, field01pixza, psi01pixza, relE01pixza, conv01pixza, niter01pixza, mallniterc01pixza, J101pixza, maxgrad01pixza, alpha01pixza, invHess01pixza] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield01pixza, 1e-3, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixza, field01pixza, psi01pixza, relE01pixza, conv01pixza, niter01pixza, mallniterc01pixza, J101pixza, maxgrad01pixza, alpha01pixza, invHess01pixza] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield01pixz, 1e-3, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixza, field01pixza, psi01pixza, relE01pixza, conv01pixza, niter01pixza, mallniterc01pixza, J101pixza, maxgrad01pixza, alpha01pixza, invHess01pixza] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield01pixz.', 1e-3, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixz, field01pixz, psi01pixz, relE01pixz, conv01pixz, niter01pixz, mallniterc01pixz, J101pixz, maxgrad01pixz, alpha01pixz, invHess01pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 10*cos(10*t), 1e-4, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixz, field01pixz, psi01pixz, relE01pixz, conv01pixz, niter01pixz, mallniterc01pixz, J101pixz, maxgrad01pixz, alpha01pixz, invHess01pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 10*cos(10*t), 1e-3, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixza, field01pixza, psi01pixza, relE01pixza, conv01pixza, niter01pixza, mallniterc01pixza, J101pixza, maxgrad01pixza, alpha01pixza, invHess01pixza] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield01pixz.', 1e-4, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixza, field01pixza, psi01pixza, relE01pixza, conv01pixza, niter01pixza, mallniterc01pixza, J101pixza, maxgrad01pixza, alpha01pixza, invHess01pixza] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield01pixz.', 1e-5, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
[allfield01pixza, field01pixza, psi01pixza, relE01pixza, conv01pixza, niter01pixza, mallniterc01pixza, J101pixza, maxgrad01pixza, alpha01pixza, invHess01pixza] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield01pixz.', 1e-4, options01pi, 0.1*pi, 1e-3*pi, 7, 7, 1e-4, 1e4);
1-J101pixza
max(abs(field01pixza))
[allfield02pixz, field02pixz, psi02pixz, relE02pixz, conv02pixz, niter02pixz, mallniterc02pixz, J102pixz, maxgrad02pixz, alpha02pixz, invHess02pixz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 5*cos(t), 1e-3, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
1-J102pixza
1-J102pixz
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, field02pixz, 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, field02pixz.', 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
options01pixz = optionsOCqn(1e-4,1e4);
options01pixz.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 200)
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-4, optionspi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-4, options01pi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-5, options01pi, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-5, options01pixz, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
options01pixz.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 500)
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-5, options01pixz, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-4, options01pixz, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
options01pixz.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 100)
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-6 6], muxz, allfield02pixz.', 1e-4, options01pixz, 0.2*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-50 50], muxz, allfield02pixz.', 1e-4, options01pixz, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixz, field02pixz, psi02pixz, relE02pixz, conv02pixz, niter02pixz, mallniterc02pixz, J102pixz, maxgrad02pixz, alpha02pixz, invHess02pixz] = OCqn([1;0], [0;1], H0, [-20 20], muxz, @(t) 5*cos(t), 1e-3, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixz, field02pixz, psi02pixz, relE02pixz, conv02pixz, niter02pixz, mallniterc02pixz, J102pixz, maxgrad02pixz, alpha02pixz, invHess02pixz] = OCqn([1;0], [0;1], H0, [-20 20], muxz, @(t) 5*cos(t), 1e-4, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-20 20], muxz, allfield02pixz, 1e-4, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-20 20], muxz, allfield02pixz.', 1e-4, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-20 20], muxz, allfield02pixz.', 1e-6, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-20 20], muxz, allfield02pixz.', 1e-5, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-100 100], muxz, allfield02pixz.', 1e-5, options01pi, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
options01pixz.f_max_alpha = @(field, direction) alpha_max_x(field, direction, 200)
[allfield02pixza, field02pixza, psi02pixza, relE02pixza, conv02pixza, niter02pixza, mallniterc02pixza, J102pixza, maxgrad02pixza, alpha02pixza, invHess02pixza] = OCqn([1;0], [0;1], H0, [-100 100], muxz, allfield02pixz.', 1e-5, options01pixz, 0.2*pi, pi/800, 7, 7, 1e-4, 1e4);
1-J102pixz
[allfieldpixz, fieldpixz, psipixz, relEpixz, convpixz, niterpixz, mallnitercpixz, J1pixz, maxgradpixz, alphapixz, invHesspixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) cos(t), 1e-4, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
1-J1pixz
figure
plot(0:pi/100:pi, fieldpixz)
[allfield05pixz, field05pixz, psi05pixz, relE05pixz, conv05pixz, niter05pixz, mallniterc05pixz, J105pixz, maxgrad05pixz, alpha05pixz, invHess05pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 2*cos(t), 1e-4, optionspi, 0.5*pi, pi/100, 7, 7, 1e-4, 1e4);
[allfield05pixz, field05pixz, psi05pixz, relE05pixz, conv05pixz, niter05pixz, mallniterc05pixz, J105pixz, maxgrad05pixz, alpha05pixz, invHess05pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 2*cos(t), 1e-4, options01pi, 0.5*pi, pi/100, 7, 7, 1e-4, 1e4);
[allfield05pixz, field05pixz, psi05pixz, relE05pixz, conv05pixz, niter05pixz, mallniterc05pixz, J105pixz, maxgrad05pixz, alpha05pixz, invHess05pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 2*cos(t), 1e-4, options01pi, 0.5*pi, pi/200, 7, 7, 1e-4, 1e4);
1-J105pixz
[allfield05pixz, field05pixz, psi05pixz, relE05pixz, conv05pixz, niter05pixz, mallniterc05pixz, J105pixz, maxgrad05pixz, alpha05pixz, invHess05pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 2*cos(t), 1e-5, options01pi, 0.5*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield05pixz, field05pixz, psi05pixz, relE05pixz, conv05pixz, niter05pixz, mallniterc05pixz, J105pixz, maxgrad05pixz, alpha05pixz, invHess05pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 2*cos(t), 1e-4, options01pi, 0.5*pi, pi/200, 7, 7, 1e-4, 1e4);
[allfield05pixz, field05pixz, psi05pixz, relE05pixz, conv05pixz, niter05pixz, mallniterc05pixz, J105pixz, maxgrad05pixz, alpha05pixz, invHess05pixz] = OCqn([1;0], [0;1], H0, [-50 50], muxz, @(t) 2*cos(t), 1e-4, options01pi, 0.5*pi, pi/400, 7, 7, 1e-4, 1e4);
1-J105pixz
max(abs(field05pixz))
save TLS2pi
figure
plot(0:pi/100:pi, fieldpixz)
figure
plot(0:pi/100:2*pi, fieldxz)
plot(0:pi/50:2*pi, fieldxz)
[allfield07pixz, field07pixz, psi07pixz, relE07pixz, conv07pixz, niter07pixz, mallniterc07pixz, J107pixz, maxgrad07pixz, alpha07pixz, invHess07pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 1.5*cos(t), 1e-3, optionspi, 2/3*pi, 2*pi/300, 7, 7, 1e-4, 1e4);
1-J107pixz
[allfield07pixz, field07pixz, psi07pixz, relE07pixz, conv07pixz, niter07pixz, mallniterc07pixz, J107pixz, maxgrad07pixz, alpha07pixz, invHess07pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 1.5*cos(t), 1e-4, optionspi, 2/3*pi, 2*pi/300, 7, 7, 1e-4, 1e4);
[allfield07pixz, field07pixz, psi07pixz, relE07pixz, conv07pixz, niter07pixz, mallniterc07pixz, J107pixz, maxgrad07pixz, alpha07pixz, invHess07pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 1.5*cos(t), 1e-3, optionspi, 2/3*pi, 2*pi/300, 7, 7, 1e-4, 1e4);
max(abs(field07pixz))
[allfield07pixz1, field07pixz1, psi07pixz1, relE07pixz1, conv07pixz1, niter07pixz1, mallniterc07pixz1, J107pixz1, maxgrad07pixz1, alpha07pixz1, invHess07pixz1] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield07pixz, 1e-4, optionspi, 2/3*pi, 2*pi/300, 7, 7, 1e-4, 1e4);
[allfield07pixz1, field07pixz1, psi07pixz1, relE07pixz1, conv07pixz1, niter07pixz1, mallniterc07pixz1, J107pixz1, maxgrad07pixz1, alpha07pixz1, invHess07pixz1] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfield07pixz.', 1e-4, optionspi, 2/3*pi, 2*pi/300, 7, 7, 1e-4, 1e4);
max(abs(field07pixz1))
1-J107pixz1
figure
plot(0:pi/300:2*pi/3, field07xz)
plot(0:pi/300:2*pi/3, field07pixz)
plot(0:2*pi/300:2*pi/3, field07pixz)
figure
plot(0:2*pi/300:2*pi/3, field07pixz1)
1-J107pixz
figure
plot(0:2*pi/300:2*pi/3, field07pixz1)
[allfieldpixz1, fieldpixz1, psipixz1, relEpixz1, convpixz1, niterpixz1, mallnitercpixz1, J1pixz1, maxgradpixz1, alphapixz1, invHesspixz1] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfieldpixz, 1e-5, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
[allfieldpixz1, fieldpixz1, psipixz1, relEpixz1, convpixz1, niterpixz1, mallnitercpixz1, J1pixz1, maxgradpixz1, alphapixz1, invHesspixz1] = OCqn([1;0], [0;1], H0, [-15 15], muxz, allfieldpixz.', 1e-5, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
figure
plot(0:pi/100:pi, fieldpixz)
hold on
plot(0:pi/100:pi, fieldpixz1)
max(abs(fieldpixz-fieldpixz1))
1-J1pixz
1-J1pixz1
clf
fieldxz1ext = [fieldxz1(1:end-1), -fieldxz1(1:end-1)]
fieldxz1ext = [fieldpixz1(1:end-1), -fieldpixz1(1:end-1)]
figure
plot(0:pi/100:2*pi, fieldpixz1ext)
clear fieldxz1ext
fieldpixz1ext = [fieldpixz1(1:end-1), -fieldpixz1(1:end-1)];
plot(0:pi/100:2*pi, fieldpixz1ext)
plot(0:pi/100:1.99*pi, fieldpixz1ext)
figure
plot(0:pi/100:1.99*pi, fieldpi1ext)
size(fieldpi1ext)
plot(0:pi/50:1.98*pi, fieldpi1ext)
ifftfieldpixz1ext = ifft(fieldpixz1ext);
figure
plot(0:10, log10(2*abs(ifftfieldpixz1ext(1:11))))
plot(0:10, (2*abs(ifftfieldpixz1ext(1:11))))
hold on
clf
plot(0:10, (2*abs(ifftfieldpixz1ext(1:11))))
hold on
plot(0:10, (2*abs(ifftfieldpi1ext(1:11))))
log(ifftfieldpi1ext(1:10)./abs(ifftfieldpi1ext(1:10)))/(-1i)/pi
log(ifftfieldpixz1ext(1:10)./abs(ifftfieldpixz1ext(1:10)))/(-1i)/pi
pi/20
angle(ifftfieldpixz1ext(1:10)./abs(ifftfieldpixz1ext(1:10)))/pi
doc fft
figure
plot(0:100, (2*abs(ifftfieldpixz1ext(1:101))))
iwfieldpixz1ext = instw(fieldpixz1ext, pi/100);
figure
plot(0:pi/100:1.99*pi, iwfieldpixz1ext)
iwfieldpixz1 = instw(fieldpixz1, pi/100);
figure
plot(0:pi/100:0.99*pi, iwfieldpixz1)
plot(0:pi/100:pi, iwfieldpixz1)
hold on
plot(0:pi/100:pi, fieldpixz1)
plot(0:pi/100:pi, fieldpixz1ext)
plot(0:pi/100:2*pi, fieldpixz1ext)
plot(0:pi/100:1.99*pi, fieldpixz1ext)
iwfieldpixz1extF = instwFourier(fieldpixz1ext, pi/100);
iwfieldpixz1extF = instwFourier(fieldpixz1ext, 2*pi);
figure
plot(0:pi/100:1.99*pi, phase)
plot(0:pi/100:1.99*pi, phi)
iw0pixz1ext = 1+2*fieldpixz1ext;
plot(0:pi/100:1.99*pi, iw0pixz1ext)
plot(0:pi/100:1.99*pi, abs(iw0pixz1ext))
save TLS2pi
%-- 18/05/2021 16:41 --%
load TLS2pi
plot(0:pi/100:1.99*pi, iw0pixz1ext)
hold on
plot(0:pi/100:1.99*pi, fieldpixz1ext)
plot(0:pi/100:pi, iwfieldpixz1)
plot(0:pi/100:1.99*pi, iwfieldpixz1ext)
figure
plot(0:pi/100:pi, iwfieldpixz1ext(1:101))
hold on
plot(0:pi/100:pi, iw0pixz1ext(1:101))
plot(0:pi/100:pi, abs(iw0pixz1ext(1:101)))
plot(0:pi/100:pi, iwfieldpixz1)
figure
plot(0:pi/50:1.98*pi, fieldpi1ext)
figure
plot(0:pi/300:2*pi/3, field07pixz)
plot(0:2*pi/300:2*pi/3, field07pixz)
hold on
plot(0:2*pi/300:2*pi/3, field07pixz1)
figure
plot(0:pi/50:2*pi, fieldxz)
max(abs(fieldxz-fieldxz(end:-1:1)))
hold on
plot(0:pi/50:2*pi, fieldxz(end:-1:1))
iwfieldxz = instw(fieldxz, pi/50);
plot(0:pi/50:2*pi, iwfieldxz)
[wrealxz, dxz, wxz, gxz, meandxz, w_maxdxz, maxdxz] = KBDMrinp(fieldxz, pi/50, 50);
[wrealxz, dxz, wxz, gxz, meandxz, w_maxdxz, maxdxz] = KBDMrinp(fieldxz, pi/50, 25);
gxz
hold on
plot(wrealxz, gxz)
wxz(1:10)
wrealxz(5)
wrealxz(5)/wrealxz(3)
figure
plot(wrealxz, meandxz)
plot(wrealxz(1:10), meandxz(1:10))
plot(wrealxz(1:5), meandxz(1:5))
plot(wrealxz(1:7), meandxz(1:7))
iwfieldxz_a = instw(fieldxz(1:end-1), pi/50);
plot(0:pi/50:2*pi, iwfieldxz)
plot(0:pi/50:1.98*pi, iwfieldxz_a)
clear iwfieldxz_a
iw0xz = 1+2*fieldxz;
plot(0:pi/50:2*pi, iw0xz)
iwfield07pixz = instw(field07pixz, 2*pi/300);
figure
plot(0:2*pi/300:2*pi/3, iwfield07pixz)
figure
whos
plot(0:pi/50:4*pi, fieldxz4)
iwfieldxz4 = instw(fieldxz4, pi/50);
hold on
plot(0:pi/50:4*pi, iwfieldxz4)
figure
plot(0:pi/50:4*pi, iwfieldxz4)
iw0xz4 = 1+2*fieldxz4;
hold on
plot(0:pi/50:4*pi, iw0xz4)
figure
plot(0:2*pi/300:2*pi/3, psi07pixz.*conj(psi07pixz))
figure
plot(0:2*pi/300:2*pi/3, psi07pixz1.*conj(psi07pixz1))
plot(0:0.01:1, fieldpixz1(1:101))
hold on
plot(1.01:0.01:2, fieldpixz1ext(102:201), '-', 'r')
plot(1.01:0.01:1.99, fieldpixz1ext(102:200), '-', 'r')
plot(1.01:0.01:1.99, fieldpixz1ext(102:200), '-')
plot(1:0.01:1.99, fieldpixz1ext(101:200), '--')
xlabel('$\frac{\tau}{\pi}$', 'interpreter', 'latex')
plot(1:0.01:2, -fieldpixz1, ':')
figure
plot(0:0.02:1, fieldpi)
hold on
plot(1:0.02:2, -fieldpi, ':')
xlabel('$\frac{\tau}{\pi}$', 'interpreter', 'latex')
figure
plot(0:1/300:1, field07pixz)
size(field07pixz)
plot(0:0.01:1, field07pixz)
hold on
plot(0:0.01:1, field07pixz1)
1-J107pixz
1-J107pixz1
xlabel('$\frac{3\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
1-J1pixz
1-J1pi
figure
plot(0:0.01:1, iwfieldpixz1ext(1:101))
hold on
plot(0:0.01:1, iw0pixz1ext(1:101))
plot(0:0.01:1, abs(iw0pixz1ext(1:101)))
xlabel('$\frac{\tau}{\pi}$', 'interpreter', 'latex')
figure
plot(0:0.01:1, fieldpixz1(1:101))
hold on
plot(0:0.01:1, iw0pixz1ext(1:101))
plot(0:0.01:1, iwfieldpixz1ext(1:101))
xlabel('$\frac{\tau}{\pi}$', 'interpreter', 'latex')
figure
plot(0:0.01:1, fieldpixz1(1:101))
muxz
figure
plot(0:pi/50:2*pi, psixz.*conj(psixz))
figure
plot(0:pi/50:2*pi, fieldxz)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
plot(0:0.02:2*pi, fieldxz)
plot(0:0.02:1, fieldxz)
plot(0:0.01:1, fieldxz)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
figure
plot(0:0.005:1, fieldxz4)
xlabel('$\frac{\tau}{4\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
figure
plot(0:pi/100:pi, psipixz1.*conj(psipixz1))
plot(0:0.01:1, psipixz1.*conj(psipixz1))
figure
plot(0:0.01:1, iwfieldxz4)
plot(0:0.005:1, iwfieldxz4)
hold on
plot(0:0.005:1, iw0xz4)
xlabel('$\frac{\tau}{4\pi}$', 'interpreter', 'latex')
figure
plot(0:0.005:1, iwfieldxz4-1)
hold on
plot(0:0.005:1, iw0xz4-1)
save TLS2pi
whos
U1genfWopt
syms n
assume(n, 'positive', 'integer')
assume(n, {'positive', 'integer'})
U1npifWopt = subs(U1genfWopt, tauf, n*pi)
U1npifWopt = simplify(subs(U1genfWopt, tauf, n*pi))
U2genfWopt = simplify(subs(U2gen, W, pi/tauf))
U2genWopt = simplify(subs(U2gen, W, pi/tauf))
U1genfWopt = simplify(subs(U1genfWopt, tau, tauf))
U2genfWopt = simplify(subs(U2genfWopt, tau, tauf))
U2npifWopt = simplify(subs(U2genfWopt, tauf, n*pi))
512/8
160/8
512/128
save TLS2pi
zeta
sigma_phi
[Psigma_phi, Dsigma_phi] = eig(sigma_phi)
sigma_phi*sigmax
sigma_phi*sigma_x
inv(Psigma_phi)
Psigma_phi = Psigma_phi/sqrt(2)
inv(Psigma_phi)
inv(Psigma_phi)-Psigma_phi'
simplify(ans)
assume(phase, 'real')
inv(Psigma_phi)-Psigma_phi'
save TLS2pi
Psigma_phi'*zeta*Psigma_phi
simplify(ans)
%-- 24/05/2021 10:05 --%
load TLS2pi
Psigma_phi
Psigma_phi'*sigma_phi*Psigma_phi
Psigma_phi'*Dsigma_phi*Psigma_phi
simplify(ans)
assume(phase, 'real')
Psigma_phi'*Dsigma_phi*Psigma_phi
sigma_phi
eig(sigma_phi)
[Psigma_phi, Dsigma_phi] = eig(sigma_phi)
Psigma_phi'*Dsigma_phi*Psigma_phi
Psigma_phi = Psigma_phi/sqrt(2)
Psigma_phi*Dsigma_phi*Psigma_phi'
Psigma_phi'*zeta*Psigma_phi
simplify(ans)
zetaD=simplify(Psigma_phi'*zeta*Psigma_phi)
exp(iW/2*D*tau)
exp(1i*W/2*D*tau)
zetaD*exp(1i*W/2*D*tau)
expm(1i*W/2*D*tau)
zetaD*expm(1i*W/2*D*tau)
expm(1i*W/2*D*tau)*zetaD*expm(1i*W/2*D*tau)
U1gen
expm(-1i*W/2*D*tau)*zetaD*expm(1i*W/2*D*tau)
M1D = expm(-1i*W/2*D*tau)*zetaD*expm(1i*W/2*D*tau)
integM1D = simplify(int(M1D, tau, 0, tau))
U1D = simplify(expm(1i*W/2*D*tau)*integM1D)
U1D = simplify(1i*W/2*expm(1i*W/2*D*tau)*integM1D)
U1
U1gen
simplify(U1gen-Psigma_phi*U1D*Psigma_phi')
M2D = simplify(expm(-1i*W/2*D*tau)*zetaD*U1D)
varphi0D = Psigma_phi'*[1;0]
integM2D = simplify(int(M2D, tau, 0, tau))
U2D = simplify(1i*W/2*expm(1i*W/2*D*tau)*integM2D)
U2gen
simplify(U2gen-Psigma_phi*U2D*Psigma_phi')
M3D = simplify(expm(-1i*W/2*D*tau)*zetaD*U2D)
integM3D = simplify(int(M3D, tau, 0, tau))
U3D = simplify(1i*W/2*expm(1i*W/2*D*tau)*integM3D)
save TLS2pi
U3DfWopt = subs(U3D, tau, tauf);
U3DfWopt = simplify(subs(U3DfWopt, W, pi/tauf))
U3DnpifWopt = simplify(subs(U3DfWopt, tauf, n*pi)
U3DnpifWopt = simplify(subs(U3DfWopt, tauf, n*pi))
assume(n, {'positive', 'integer'})
U3DnpifWopt = simplify(subs(U3DfWopt, tauf, n*pi))
save TLS2pi
%-- 26/05/2021 17:58 --%
load TLS2pi
whos
muxz1 = [0.5 1; 1 -0.5]
muxz05 = [0.5 1; 1 -0.5]
clear muxz1
[allfieldpixz05, fieldpixz05, psipixz05, relEpixz05, convpixz05, niterpixz05, mallnitercpixz05, J1pixz05, maxgradpixz05, alphapixz05, invHesspixz05] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, @(t) cos(t), 1e-3, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
1-J1pixz05
figure
plot(0:pi/100:pi, fieldpixz05)
1-J1pixz1
1-J1pixz
[allfieldpixz05a, fieldpixz05a, psipixz05a, relEpixz05a, convpixz05a, niterpixz05a, mallnitercpixz05a, J1pixz05a, maxgradpixz05a, alphapixz05a, invHesspixz05a] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, allfieldpixz.', 1e-4, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
[allfieldpixz05a, fieldpixz05a, psipixz05a, relEpixz05a, convpixz05a, niterpixz05a, mallnitercpixz05a, J1pixz05a, maxgradpixz05a, alphapixz05a, invHesspixz05a] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, allfieldpixz05.', 1e-4, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
[allfieldpixz05a, fieldpixz05a, psipixz05a, relEpixz05a, convpixz05a, niterpixz05a, mallnitercpixz05a, J1pixz05a, maxgradpixz05a, alphapixz05a, invHesspixz05a] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, allfieldpixz05.', 1e-5, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
1-J1pixz05a
[allfieldpixz05a, fieldpixz05a, psipixz05a, relEpixz05a, convpixz05a, niterpixz05a, mallnitercpixz05a, J1pixz05a, maxgradpixz05a, alphapixz05a, invHesspixz05a] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, @(t) cos(t), 1e-4, optionspi, pi, pi/100, 7, 7, 1e-4, 1e4);
1-J1pixz05a
figure
plot(0:pi/100:pi, fieldpixz05a)
hold on
plot(0:pi/100:pi, fieldpixz05a)
fieldpixz05aext = [fieldpixz05a(1:end-1), -fieldpixz05a(1:end-1)];
figure
plot(0:pi/100:1.99*pi, fieldpixz05aext)
ifftfieldpixz05aext = ifft(fieldpixz05aext);
figure
plot(0:100, (2*abs(ifftfieldpixz05ext(1:101))))
plot(0:100, (2*abs(ifftfieldpixz05aext(1:101))))
hold on
plot(0:100, (2*abs(ifftfieldpiext(1:101))))
plot(0:100, (2*abs(ifftfieldpi1ext(1:101))))
plot(0:50, (2*abs(ifftfieldpi1ext(1:51))))
hold on
plot(0:pi/50:pi, fieldpi1)
plot(0:pi/50:1.98*pi, fieldpi1ext)
iwfieldxz05aext = instw(fieldxz05aext, pi/100);
iwfieldxz05aext = instw(fieldpixz05aext, pi/100);
hold on
plot(0:pi/100:1.99*pi, iwfieldxz05aext)
[allfield05pixz05, field05pixz05, psi05pixz05, relE05pixz05, conv05pixz05, niter05pixz05, mallniterc05pixz05, J105pixz05, maxgrad05pixz05, alpha05pixz05, invHess05pixz05] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, @(t) cos(t), 1e-3, optionspi, pi/2, pi/100, 7, 7, 1e-4, 1e4);
1-J105pixz05
figure
plot(0:pi/100:pi/2, field05pixz05)
[allfield05pixz05a, field05pixz05a, psi05pixz05a, relE05pixz05a, conv05pixz05a, niter05pixz05a, mallniterc05pixz05a, J105pixz05a, maxgrad05pixz05a, alpha05pixz05a, invHess05pixz05a] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, allfield05pixz05.', 1e-4, optionspi, pi/2, pi/100, 7, 7, 1e-4, 1e4);
1-J105pixz05a
figure
plot(0:pi/100:pi/2, field05pixz05a)
hold on
plot(0:pi/100:pi/2, field05pixz05a)
figure
plot(0:pi/100:pi, fieldpixz1)
[allfield03pixz05, field03pixz05, psi03pixz05, relE03pixz05, conv03pixz05, niter03pixz05, mallniterc03pixz05, J103pixz05, maxgrad03pixz05, alpha03pixz05, invHess03pixz05] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, @(t) 3*cos(t), 1e-3, optionspi, pi/3, pi/300, 7, 7, 1e-4, 1e4);
1-J103pixz05
[allfield03pixz05, field03pixz05, psi03pixz05, relE03pixz05, conv03pixz05, niter03pixz05, mallniterc03pixz05, J103pixz05, maxgrad03pixz05, alpha03pixz05, invHess03pixz05] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, @(t) 3*cos(t), 1e-4, optionspi, pi/3, pi/300, 7, 7, 1e-4, 1e4);
optionspi
[allfield03pixz05, field03pixz05, psi03pixz05, relE03pixz05, conv03pixz05, niter03pixz05, mallniterc03pixz05, J103pixz05, maxgrad03pixz05, alpha03pixz05, invHess03pixz05] = OCqn([1;0], [0;1], H0, [-15 15], muxz05, @(t) 3*cos(t), 1e-4, options01pi, pi/3, pi/300, 7, 7, 1e-4, 1e4);
1-J103pixz05
figure
plot(0:pi/300:pi/3, field03pixz05)
pi/20
[allfield03pixz05, field03pixz05, psi03pixz05, relE03pixz05, conv03pixz05, niter03pixz05, mallniterc03pixz05, J103pixz05, maxgrad03pixz05, alpha03pixz05, invHess03pixz05] = OCqn([1;0], [0;1], H0, [-60 60], muxz05, @(t) 3*cos(t), 1e-4, options01pi, pi/3, pi/300, 7, 9, 1e-4, 1e4);
figure
plot(0:pi/300:pi/3, field03pixz05)
1-J103pixz05
figure
plot(0:pi/100:pi, fieldpixz1)
hold on
plot(0:pi/100:pi, iwfieldpixz1)
plot(0:pi/100:pi, iwfieldpixz1ext(1:101))
save TLS2pi
figure
plot(0:0.01:1, fieldpixz05)
hold on
plot(1:0.01:2, -fieldpixz05)
xlabel('$\frac{\tau}{\pi}$', 'interpreter', 'latex')
figure
plot(0:0.01:1, fieldpixz05a)
figure
plot(0:0.01:1, field05pixz05a)
plot(0:0.02:1, field05pixz05a)
hold on
plot(1:0.02:2, -field05pixz05a)
xlabel('$\frac{2\tau}{\pi}$', 'interpreter', 'latex')
whos
figure
plot(0:0.01:2, field01pi1ext)
plot(0:0.01:1.99, field01pi1ext)
figure
plot(0:100, (2*abs(ifftfieldpi1ext(1:101))))
plot(0:99, (2*abs(ifftfieldpi1ext(1:100))))
plot(0:20, (2*abs(ifftfieldpi1ext(1:21)))
plot(0:20, (2*abs(ifftfieldpi1ext(1:21)))))
plot(0:20, (2*abs(ifftfieldpi1ext(1:21))))
plot(0:20, (2*abs(ifftfieldpixz1ext(1:21))))
save TLS2pi
[allfield8pixz, field8pixz, psi8pixz, relE8pixz, conv8pixz, niter8pixz, mallniterc8pixz, J18pixz, maxgrad8pixz, alpha8pixz, invHess8pixz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 1/8*cos(t), 1e-4, options, 8*pi, pi/25, 7, 7, 1e-4, 1e4);
[allfield8pixz, field8pixz, psi8pixz, relE8pixz, conv8pixz, niter8pixz, mallniterc8pixz, J18pixz, maxgrad8pixz, alpha8pixz, invHess8pixz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 1/8*cos(t), 1e-4, options, 8*pi, pi/50, 7, 7, 1e-4, 1e4);
[allfield8pixz, field8pixz, psi8pixz, relE8pixz, conv8pixz, niter8pixz, mallniterc8pixz, J18pixz, maxgrad8pixz, alpha8pixz, invHess8pixz] = OCqn([1;0], [0;1], H0, [-15 15], muxz, @(t) 1/8*cos(t), 1e-4, options, 8*pi, pi/50, 7, 7, 1e-4, 1e4);
[allfield8pixz, field8pixz, psi8pixz, relE8pixz, conv8pixz, niter8pixz, mallniterc8pixz, J18pixz, maxgrad8pixz, alpha8pixz, invHess8pixz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 1/8*cos(t), 1e-3, options, 8*pi, pi/50, 7, 7, 1e-4, 1e4);
[allfield8pixz, field8pixz, psi8pixz, relE8pixz, conv8pixz, niter8pixz, mallniterc8pixz, J18pixz, maxgrad8pixz, alpha8pixz, invHess8pixz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 1/8*cos(t), 1e-3, options, 8*pi, pi/100, 7, 7, 1e-4, 1e4);
1-J18pixz
figure
plot(0:pi/100:8*pi, field8pixz)
[allfield16pixz, field16pixz, psi16pixz, relE16pixz, conv16pixz, niter16pixz, mallniterc16pixz, J116pixz, maxgrad16pixz, alpha16pixz, invHess16pixz] = OCqn([1;0], [0;1], H0, [-6 6], muxz, @(t) 1/16*cos(t), 1e-3, options, 16*pi, pi/50, 7, 7, 1e-4, 1e4);
1-J116pixz
figure
plot(0:pi/50:16*pi, field16pixz)
1/16
figure
plot(0:pi/50:4*pi, field4pixz)
plot(0:pi/50:4*pi, field4xz)
whos
plot(0:pi/50:4*pi, fieldxz4)
save TLS2pi
whos
clear invHess8pixz invHess16pixz
save TLS2pi
%-- 14/06/2021 11:24 --%
doc ifft
whos
load TLS2pi
whos
Afield_w1
figure
plot(phases/(pi)-0.5, allJmax)
figure
plot(phases/(pi)-0.5, allJmax_maxf)
1-max(allJmax_maxf)
figure
plot(phases/(pi)-0.5, allJmax_Afield_w1)
plot(phases/(pi)-0.5, allJmax_Aw1)
ifftfield(1:10)
abs(ifftfield(1:10))
[abs(ifftfield(1),  2*abs(ifftfield(2:10))]
[abs(ifftfield(1),  2*abs(ifftfield(2:10)))]
[abs(ifftfield(1)),  2*abs(ifftfield(2:10))]
-arg(1:10)/1i
-angle(ifftfield(1:10))
-angle(ifftfield(1:10))/pi
cfield  = [ifftfield(1), 2*ifftfield(2:50)];
figure
plot(0:9, cfield(1:10))
cfield  = [abs(ifftfield(1)), 2*abs(ifftfield(2:50))];
plot(0:9, cfield(1:10))
plot(0:9, log10(cfield(1:10)))
plot(0:10, log10(cfield(1:11)))
xlabel('$n$', 'interpreter', 'latex')
ylabel('$c_n$', 'interpreter', 'latex')
ylabel('$\log_{10}(c_n)$', 'interpreter', 'latex')
allJmax_Aw1fine(615)
331*2*pi/1e3
331*2/1e3 - 0.5
340*2/1e3 - 0.5
332*2/1e3 - 0.5
figure
plot((331:2e-3:333)*2/1e3 - 0.5, allJmax_Aw1fine)
allJmax_Aw1fine(10)
1-allJmax_Aw1fine(10)
1-max(allJmax_Aw1fine)
figure
plot(0:5e-3:1, field4pi)
hold on
plot(0:5e-3:1, field4pi1)
1-J14pi
1-J14pi1
cfield4pi = cfield  = [ifftfield4pi(1), 2*ifftfield4pi(2:100)];
cfield4pi = [ifftfield4pi(1), 2*ifftfield4pi(2:100)];
cfield4pi
cfield4pi = abs([ifftfield4pi(1), 2*ifftfield4pi(2:100)]);
cfield4pi(1:11)
cfield4pi(1:21)
figure
plot(0:0.5:10, cfield4pi(1:21))
plot(0:0.5:10, log10(cfield4pi(1:21)))
figure
plot(0:0.5:50, log10(cfield4pi))
plot(0:0.5:49.5, log10(cfield4pi))
-angle(ifftfield(1:21))/pi
-angle(ifftfield4pi(1:21))/pi
-angle(ifftfield4pi)/pi
max(abs(field4pi-field4pi(end:-1:1)))
max(abs(field4pi(1:100) + field4pi(101:200)))
size(field4pi)
max(abs(field4pi(1:150) + field4pi(51:200)))
hold on
plot(0:5e-3:0.75, -field4pi1(51:200))
plot(0:5e-3:0.745, -field4pi(51:200))
plot(0:5e-3:0.75, -field4pi(50:200))
plot(0:5e-3:0.745, -field4pi(49:199))
plot(0:5e-3:0.745, -field4pi(50:199))
2-1.86/pi
pi-ans
2-1.86
max(abs(fieldxz-fieldxz(end:-1:1)))
max(abs(fieldxz-fieldxz(end:-1:1)))/max(abs(fieldxz))
max(abs(field01pi1))
max(abs(field02pi2))
max(abs(field02pi2))/5
max(abs(field02pi1))/5
max(abs(field01pi1))
max(abs(fieldxzpi))
max(abs(fieldpixz))
max(abs(fieldxz))
max(abs(fieldxz))/0.5
figure
size(field4pi)
plot(0:5e-3:1, field4pi1)
plot(0:5e-3:1, field4pi)
xlabel('$\frac{\tau}{4\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
whos
figure
plot(0:0.025:1, field01pi1)
size(field01pi1)
plot(0:0.025:1, field02pi2)
size(field02pi2)
plot(0:0.0125:1, field02pi2)
hold on
plot(1:0.0125:2, -field02pi2)
xlabel('$\frac{5\tau}{\pi}$', 'interpreter', 'latex')
figure
plot(0:0.01:1, field01pi1)
hold on
plot(1:0.01:2, -field01pi1)
xlabel('$\frac{10\tau}{\pi}$', 'interpreter', 'latex')
1-J1field07xz
1-J1field07pixz
1-J107pixz
1-J107pixz1
1-J107pixz2
max(abs(field07pixz - field07pixz(end:-1:1))
max(abs(field07pixz - field07pixz(end:-1:1)))
max(abs(field07pixz1 - field07pixz1(end:-1:1)))
max(abs(field07pixz - field07pixz(end:-1:1)))./max(abs(field07pixz))
max(abs(field07pixz1 - field07pixz1(end:-1:1)))/max(abs(field07pixz1))
max(abs(fieldpixz05))
max(abs(field05pixz05)
max(abs(field05pixz05))
max(abs(field05pixz05))/2
%-- 28/06/2021 16:08 --%
load TLS2pi
whos
%-- 29/06/2021 9:59 --%
load TLS2pi
whos
figure
plot(0:pi/50:2*pi, field)
figure
plot(0:pi/100:pi, fieldpixz1)
save epsilon_data field fieldpixz1
load('epsilon_data.mat')
muxz_exp = [1/sqrt(2) 1; 1 -1/sqrt(2)]
[allfieldexpxz, fieldexpxz, psiexpxz, relEexpxz, convexpxz, niterexpxz, mallnitercexpxz, J1expxz, maxgradexpxz, alphaexpxz, invHessexpxz] = OCqn([1;0], [0;1], H0, [-6 6], muxz_exp, @(t) 0.5*cos(t), 1e-4, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
figure
plot(0:pi/50:2*pi, fieldexpxz)
hold on
plot(0:pi/50:2*pi, fieldxz1)
plot(0:pi/50:2*pi, fieldxz)
1-J1expxz
[allfieldexpxzlim, fieldexpxzlim, psiexpxzlim, relEexpxzlim, convexpxzlim, niterexpxzlim, mallnitercexpxzlim, J1expxzlim, maxgradexpxzlim, alphaexpxzlim, invHessexpxzlim] = OCqn([1;0], [0;1], H0, [-6 6], muxz_exp, @(t) 0.5*cos(t), 1e-2, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
figure
plot(0:pi/50:2*pi, fieldexpxzlim)
1-J1expxzlim
plot(0:pi/50:2*pi, fieldexpxzlim)
[allfieldexpxzlim, fieldexpxzlim, psiexpxzlim, relEexpxzlim, convexpxzlim, niterexpxzlim, mallnitercexpxzlim, J1expxzlim, maxgradexpxzlim, alphaexpxzlim, invHessexpxzlim] = OCqn([1;0], [0;1], H0, [-6 6], muxz_exp, @(t) 0.5*cos(t), 1e-1, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1expxzlim
plot(0:pi/50:2*pi, fieldexpxzlim)
[allfieldexpxzlim1, fieldexpxzlim1, psiexpxzlim1, relEexpxzlim1, convexpxzlim1, niterexpxzlim1, mallnitercexpxzlim1, J1expxzlim1, maxgradexpxzlim1, alphaexpxzlim1, invHessexpxzlim1] = OCqn([1;0], [0;1], H0, [-6 6], muxz_exp, @(t) 0.5*cos(t), 1, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1expxzlim1
[allfieldexpxzlim1, fieldexpxzlim1, psiexpxzlim1, relEexpxzlim1, convexpxzlim1, niterexpxzlim1, mallnitercexpxzlim1, J1expxzlim1, maxgradexpxzlim1, alphaexpxzlim1, invHessexpxzlim1] = OCqn([1;0], [0;1], H0, [-6 6], muxz_exp, @(t) 0.5*cos(t), 0.2, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1expxzlim1
1-J1expxzlim
figure
plot(0:pi/50:2*pi, fieldexpxzlim1)
[allfieldexpxzlim2, fieldexpxzlim2, psiexpxzlim2, relEexpxzlim2, convexpxzlim2, niterexpxzlim2, mallnitercexpxzlim2, J1expxzlim2, maxgradexpxzlim2, alphaexpxzlim2, invHessexpxzlim2] = OCqn([1;0], [0;1], H0, [-6 6], muxz_exp, @(t) 0.5*cos(t), 0.3, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-4, 1e4);
1-J1expxzlim2
plot(0:pi/50:2*pi, fieldexpxzlim2)
save TLS2pi
figure
plot(0:0.01:1, fieldexpxz)
hold on
plot(0:0.01:1, fieldexpxzlim)
plot(0:0.01:1, fieldexpxzlim1)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
1-J1expxz
1-J1expxzlim
1-J1expxzlim1
save TLS2pi
save epsilon_exp fieldexpxz fieldexpxzlim fieldexpxzlim1
%-- 08/07/2021 18:22 --%
load TLS2pi
cftool
nu_cut = 26.64/2.5
stiffness = 0.1155*2.5
filterE = @(nu) 0.5*(1-tanh(stiffness*(nu-nu_cut)))
figure
plot(0:0.5:50, filterE(0:0.5:50))
filterEcut = @(nu) 0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20))
figure
plot(0:0.5:50, filterEcut(0:0.5:50))
filterEcut = @(nu) 0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1))
plot(0:0.5:50, filterEcut(0:0.5:50))
[Hoperations_exp, fcouplingOp_exp] = Hmats2Hops(H0, mu)
H0
muxz_exp
[Hoperations_exp, fcouplingOp_exp] = Hmats2Hops(H0, muxz_exp)
[fieldtf, fieldnuf, psif, relEf, convf, niterf, mallnitercf, Jtermsf, maxgradf, alphaf, invHessf] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], fguess, filterEcut, options, 2*pi, 2*pi/1e2, 7, 7, 1e-6);
ftarget = get_fchi_state([0;1])
figure
plot(0:pi/50:2*pi, cos(0:pi/50:2*pi))
dctIcos2pi = dctI(0.5*cos(0:pi/50:2*pi))*2*sqrt(pi/100);
figure
plot(0:0.5:50, dctIcos2pi)
dctIcos2pi(abs(dctIcos2pi)<1e-8)= 0;
plot(0:0.5:50, dctIcos2pi)
[fieldtf, fieldnuf, psif, relEf, convf, niterf, mallnitercf, Jtermsf, maxgradf, alphaf, invHessf] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), filterEcut, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-6);
[fieldtf, fieldnuf, psif, relEf, convf, niterf, mallnitercf, Jtermsf, maxgradf, alphaf, invHessf] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41).', filterEcut, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
[fieldtf, fieldnuf, psif, relEf, convf, niterf, mallnitercf, Jtermsf, maxgradf, alphaf, invHessf] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), filterEcut, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
filterEcut = @(nu) 0.1*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1))
[fieldtf, fieldnuf, psif, relEf, convf, niterf, mallnitercf, Jtermsf, maxgradf, alphaf, invHessf] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), filterEcut, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
filterEcut = @(nu) 10*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1))
[fieldtf, fieldnuf, psif, relEf, convf, niterf, mallnitercf, Jtermsf, maxgradf, alphaf, invHessf] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), filterEcut, optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
figure
plot(0:pi/50:2*pi, fieldtf)
1-J1f
Jtermsf
figure
plot(0:0.5:20, fieldnuf)
length(fieldnuf)
plot(0:0.5:31, fieldnuf(1:31))
plot(0:0.5:30, fieldnuf(1:61))
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
Jtermsf1
figure
plot(0:0.5:20, fieldnuf(1:41))
plot(0:0.5:20, fieldnuf1(1:41))
figure
plot(0:pi/50:2*pi, fieldtf1)
figure
plot(0:pi/50:2*pi, psif.*conj(psif))
plot(0:pi/50:2*pi, psif.*conj(psif1))
plot(0:pi/50:2*pi, psif1.*conj(psif1))
Jtermsf
Jtermsf1
1-Jtermsf1.Jmax
1-Jtermsf.Jmax
[fieldtf2, fieldnuf2, psif2, relEf2, convf2, niterf2, mallnitercf2, Jtermsf2, maxgradf2, alphaf2, invHessf2] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 50*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
1-Jtermsf2.Jmax
figure
plot(0:pi/50:2*pi, fieldtf2)
[fieldtf2, fieldnuf2, psif2, relEf2, convf2, niterf2, mallnitercf2, Jtermsf2, maxgradf2, alphaf2, invHessf2] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 25*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
1-Jtermsf2.Jmax
figure
plot(0:pi/50:2*pi, fieldtf2)
[fieldtf2, fieldnuf2, psif2, relEf2, convf2, niterf2, mallnitercf2, Jtermsf2, maxgradf2, alphaf2, invHessf2] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 40*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
1-Jtermsf2.Jmax
plot(0:pi/50:2*pi, fieldtf2)
plot(0:pi/50:2*pi, fieldtf1)
plot(0:0.01:1, fieldtf1)
hold on
plot(0:0.01:1, fieldtf)
clf
plot(0:0.01:1, fieldtf)
plot(0:0.01:1, fieldtf1)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
save epsilon_exp_limf fieldtf1
fftfieldtf1 = fft(fieldtf1);
figure
plot(0:100, fftfieldtf1)
plot(0:100, 2*abs(fftfieldtf1))
iwfieldf1 = instw(fieldf1, pi/50);
iwfieldf1 = instw(fieldtf1, pi/50);
figure
plot(0:pi/50:2*pi, iwfieldf1)
iwcosfieldf1 = instwcos(fieldtf1, 2*pi);
hold on
plot(0:pi/50:2*pi, iwcosfieldf1)
hold on
plot(0:pi/50:2*pi, 1+2*fieldtf1)
fftfieldtf1 = fft(fieldtf1(1:end-1));
plot(0:99, 2*abs(fftfieldtf1))
[fieldtf2, fieldnuf2, psif2, relEf2, convf2, niterf2, mallnitercf2, Jtermsf2, maxgradf2, alphaf2, invHessf2] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 5*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
1-Jtermsf2.Jmax
figure
plot(0:0.01:1, fieldtf2)
[fieldtf2, fieldnuf2, psif2, relEf2, convf2, niterf2, mallnitercf2, Jtermsf2, maxgradf2, alphaf2, invHessf2] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 2*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
1-Jtermsf2.Jmax
plot(0:0.01:1, fieldtf2)
save TLS2pi
%-- 19/07/2021 17:30 --%
load TLS2pi
addpath 'C:\Users\idosc\Dropbox\MATLAB\OCT\QuantumComputation\TLSshort'
load TLS2pi
whos
simplify(U1D(:,2)-[conj(U1D(2,1)); -conj(U1D(2,1))])
simplify(ans)
simplify(U1D(:,2)-[conj(U1D(2,1)); -conj(U1D(1,1))])
simplify(ans)
U1D
simplify(U1D(1,1)+ conj(U1D(2,2)))
assume(W, {'real'})
assume(tau, {'real'})
simplify(U1D(1,1)+ conj(U1D(2,2)))
assume(phase, {'real'})
simplify(U1D(1,1)+ conj(U1D(2,2)))
simplify(U1D(1,1)- conj(U1D(2,2)))
simplify(U1D(2,1)+ conj(U1D(1,2)))
simplify(U2D(2,1)+ conj(U2D(1,2)))
simplify(U3D(2,1)+ conj(U3D(1,2)))
simplify(U3D(1,1)- conj(U3D(2,2)))
%-- 18/08/2021 12:38 --%
(-100+1000*log(1.25))/1.6
1/0.8
%-- 23/08/2021 18:10 --%
syms a b A B x xtag
assume(a, b, A, B, x, xtag, {'real'})
assume(a, {'real'})
doc assume
assume(b, {'real'})
assume(A, {'real'})
assume(A, {'clear'})
doc int
int(exp(1i*a*xtag)*cos(b*(xtag+phi))*xtag^n, xtag, 0, x)
syms phi n
assume(phi, {'real'})
assume(n, {'integer'})
int
int(exp(1i*a*xtag)*cos(b*(xtag+phi))*xtag^n, xtag, 0, x)
int(cos(b*(xtag+phi))*xtag^n, xtag, 0, x)
int(xtag^n, xtag, 0, x)
%-- 17/11/2021 13:30 --%
load TLS2pi
whos
U1gen
U1genfWopt
U2gen
U2genfWopt
U2npifWopt
U2D
U2D2pi
U2DW2pi