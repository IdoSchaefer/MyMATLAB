%-- 17/11/2020 11:07 --%
load HTLSs_harmonic
Hop10nl = @(v, field) H3harcus10nl*v - field(1)*(H3LSe_1hus*v) - field(2)*(H3LSe_2hus*v)
psi4cPWC40 = SchrPWCcheb(Hop10nl, u0_3, fieldtg4c(:, 2:2:end-1), [-15 15], 6, 40, 20);
[Jmax_overlap4cPWC40, phases4cPWC40] = Uoverlap_gate(psi4cPWC40(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlap4cPWC40
psi4cPWC20 = SchrPWCcheb(Hop10nl, u0_3, fieldtg4c(:, 3:4:end-1), [-15 15], 6, 20, 40);
[Jmax_overlap4cPWC20, phases4cPWC20] = Uoverlap_gate(psi4cPWC20(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlap4cPWC20
doc size
figure
plot(1:4, 2:2:8)
plot((1:4).', 2:2:8)
plot((1:4), (2:2:8).')
plot((1:4), [(2:2:8);)
plot((1:4), [(2:2:8);3:3:12])
plot((1:4).', [(2:2:8);3:3:12])
plot((1:4), [(2:2:8);3:3:12].')
plot_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1))
hold on
plot(0:0.075:6, fieldtg4c)
figure
plot_steps(0:6/40:6, fieldtg4c(:, 2:2:end-1))
hold on
plot(0:0.075:6, fieldtg4c)
erf(inf)
erf(-inf)
erf(1/sqrt(2))
erf(2/sqrt(2))
softrectfun_erf(1, 0, 2, 1/sqrt(2))
softrectfun_erf(1, 0, 2, 1/10)
x=0:0.01:10;
figure
plot(x, softrectfun_erf(1, 0, 2, 1/10)
plot(x, softrectfun_erf(x, 0, 2, 1/10)
plot(x, softrectfun_erf(x, 3, 7, 1/sqrt(2))
plot(x, softrectfun_erf(x, 3, 7, 1/sqrt(2)))
plot(x, softrectfun_erf(x, 3, 7, 4))
plot(x, softrectfun_erf(x, 3, 7, 2))
plot(x, softrectfun_erf(x, 4, 6, 2))
plot(x, softrectfun_erf(x, 4, 6, 1/10))
plot(x, softrectfun_erf(x, 4, 6, 1/3))
plot(x, softrectfun_erf(x, 4, 6, 2))
plot(x, softrectfun_erf(x, 4, 6, 1))
plot(x, softrectfun_erf(x, 4, 6, 2))
tsoft = 0:0.01:6;
f_smooth = soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), tsoft, @(x, a, b) softrectfun_erf(x, a, b, 1/3));
figure
plot(tsoft, f_
plot(tsoft, f_)3sm6
plot(tsoft, f_smooth)
plot_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1))
hold on
plot(tsoft, f_smooth)
f_smooth = soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), tsoft, @(x, a, b) softrectfun_erf(x, a, b, 1/10));
plot(tsoft, f_smooth)
6/25
f_smooth = soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), tsoft, @(x, a, b) softrectfun_erf(x, a, b, 6/25));
figure
plot(tsoft, f_smooth)
plot_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1))
hold on
plot(tsoft, f_smooth)
f_smooth1 = soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), tsoft, @(x, a, b) softrectfun_erf(x, a, b, 6/50));
plot(tsoft, f_smooth1)
6/(2*pi)*50
(6/(2*pi)*50)/2
(2*pi)/50
(2*pi)/50*2
figure
plot(x, softrectfun_erf(x(301:701), 4, 6, 2*pi/50))
plot(x(301:701), softrectfun_erf(x(301:701), 4, 6, 2*pi/50))
plot(x(301:701), softrectfun_erf(x(301:701), 5-pi/50, 5+pi/50, 2*pi/50))
fEsmooth = @(t) soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), t, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/25));
[psi4c_smooth, mniter, matvecs, max_errors, history] = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEsmooth, [], [-15 15], u0_3, 0:0.075:6, 80, 9, 7, 1e-6, 10, 16);
max_errors
[Jmax_overlap4c_smooth, phases4c_smooth] = Uoverlap_gate(psi4c_smooth(:, end), target3LSu(:), [1 4 7])
[psi4c_smooth1, mniter, matvecs, max_errors, history] = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEsmooth1, [], [-15 15], u0_3, 0:0.075:6, 80, 9, 7, 1e-6, 10, 16);
fEsmooth1 = @(t) soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), t, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/25));
[psi4c_smooth1, mniter, matvecs, max_errors, history] = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEsmooth1, [], [-15 15], u0_3, 0:0.075:6, 80, 9, 7, 1e-6, 10, 16);
[Jmax_overlap4c_smooth1, phases4c_smooth1] = Uoverlap_gate(psi4c_smooth1(:, end), target3LSu(:), [1 4 7])
fEsmooth1 = @(t) soft_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1), t, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/50));
[psi4c_smooth1, mniter, matvecs, max_errors, history] = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEsmooth1, [], [-15 15], u0_3, 0:0.075:6, 80, 9, 7, 1e-6, 10, 16);
[Jmax_overlap4c_smooth1, phases4c_smooth1] = Uoverlap_gate(psi4c_smooth1(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlap4c_smooth1
whos
save HTLSs_harmonic
test_step_smoothing
1 - allJmax_overlap4c_smooth
1 - Jmax_overlap4c_smooth
test_step_smoothing
1 - allJmax_overlap4c_smooth
plo
figure
plot((1:Nsigma)*Dsigma, log10(allJmax_overlap4c_smooth))
plot((1:Nsigma)*Dsigma, log10(1 - allJmax_overlap4c_smooth))
xlabel('$\sigma \text{(ns)}$', 'interpreter', 'latex')
xlabel('$\sigma \textrm{(ns)}$', 'interpreter', 'latex')
ylabel('$\log_{10}(infidelity)$',  'interpreter', 'latex')
hold on
plot([0 2], log10(1-Jmax_overlap4cPWC20)*[1, 1])
6/(2*pi)*50
6/(2*pi)*50/20
1 - allJmax_overlap4c_smooth
test_step_smoothing
1 - allJmax_overlap4c_smooth
1 - Jmax_overlap4c
1 - Jmax_overlapg4c
1-Jmax_overlap4cPWC20
test_step_smoothing
figure
plot_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1))
hold on
vandermonde
v
v(1:4)-v(1:4).'
[finterp, coefs] = softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), t, fxmid, xresult, @(x, a, b) softrectfun_erf(x, a, b, 2));
[finterp, coefs] = softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), 0:0.01:6, fxmid, xresult, @(x, a, b) softrectfun_erf(x, a, b, 2));
[finterp, coefs] = softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), 0:0.01:6, @(x, a, b) softrectfun_erf(x, a, b, 2));
[finterp, coefs] = softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), 0:0.01:6, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/25));
figure
plot(0:0.075:6, fieldtg4c)
hold on
plot(0:0.01:6, finterp)
clear finterp coefs
[fieldtg4cInt2, coefs] = softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), 0:0.075:6, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/25));
clf
plot(0:0.075:6, fieldtg4c)
plot(0:0.075:6, fieldtg4cInt2)
plot(0:0.075:6, fieldtg4c)
hold on
plot(0:0.075:6, fieldtg4cInt2)
max(abs(fieldtg4c - fieldtg4cInt2), [], 2)
max(abs(fieldtg4c - fieldtg4cInt2)./abs(fieldtg4c), [], 2)
max(abs(fieldtg4c - fieldtg4cInt2)./max(abs(fieldtg4c), [], 2), [], 2)
fEsmoothInt2 = @(t) soft_steps(0:6/20:6, coefs.', t, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/25));
psi4cInt2 = SemiGlobalHparams(Hoperations10nl.psi, Hoperations10nl.diff_psi, 2, fEsmoothInt2, [], [-15 15], u0_3, [0, 6], 80, 9, 7,...
1e-6, 10, 16, [], false);
[Jmax_overlap4cInt2, phases4cInt2] = Uoverlap_gate(psi4cInt2(:, end), target3LSu(:), [1 4 7])
1-Jmax_overlap4cInt2
1-Jmax_overlap4c
1-Jmax_overlapg4c
figure
save HTLSs_harmonic
plot(0:6/20:6, fieldtg4c(:, 3:4:end-1))
plot_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1))
hold on
plot_steps(0:6/20:6, coefs)
clf
plot(6/40:6/20:6-6/40, fieldtg4c(:, 3:4:end-1))
hold on
plot(6/40:6/20:6-6/40, coefs)
test_soft_int
hold on
plot([0 2], log10(1-Jmax_overlap4cPWC20)*[1, 1])
plot([0 2], log10(1-Jmax_overlapg4c)*[1, 1])
test_soft_int
hold on
plot([0 2], log10(1-Jmax_overlapg4c)*[1, 1])
plot([0 5], log10(1-Jmax_overlapg4c)*[1, 1])
plot([0 5], log10(1-Jmax_overlap4cPWC20)*[1, 1])
xlabel('$\sigma \textrm{(ns)}$', 'interpreter', 'latex')
ylabel('$\log_{10}(infidelity)$',  'interpreter', 'latex')
0.25/0.03
ans/(2*pi)
1/ans
2*pi/25
[fieldtg4cInt05, coefs05] = softrect_interp(0:6/20:6, fieldtg4c(:, 3:4:end-1), 0:0.075:6, @(x, a, b) softrectfun_erf(x, a, b, 2*pi/100));
figure
plot(0:0.075:6, fieldtg4c)
hold on
plot(0:0.075:6, fieldtg4cInt05)
plot_steps(0:6/20:6, fieldtg4c(:, 3:4:end-1))
save HTLSs_harmonic
figure
plot(x, softrectfun_erf(x, 5-pi/50, 5+pi/50, 2*pi/50))
x1 = -5:0.01:5;
plot(x1, softrectfun_erf(x, -pi/50, pi/50, 2*pi/50))
plot(x1, softrectfun_erf(x1, -pi/50, pi/50, 2*pi/50))
plot(x1(201:801), softrectfun_erf(x1(201:801), -pi/50, pi/50, 2*pi/50))
plot(x1(201:801), softrectfun_erf(x1(201:801), -1, 1, 2))
plot(x1, softrectfun_erf(x1, -1, 1, 2))
hold on
plot(x1, softrectfun_erf(x1, -1, 1, 0.1))
plot(x1, softrectfun_erf(x1, -1, 1, 0.8))
plot(x1, softrectfun_erf(x1, -1, 1, 0.5))
figure
plot_steps(0:6/20:6, fieldtg4c(1, 3:4:end-1))
hold on
plot(tsoft, f_smooth1)
plot(tsoft, fEsmooth1(tsoft))
%-- 22/11/2020 19:37 --%
doc mod
ceil(4.1)
kron_index = kron_index([2 4], [4 2 3])
inv_kron_index([2 4], 31)
kron_index = Mkron_index([2 4], [4 2 3])
kron_index = kron_index([2 4], [4 2 3; 6 1 4])
clear kron_index
kron_index([2 4], [4 2 3; 6 1 4])
inv(kron_index([2 4], ans)
inv_kron_index([2 4], ans)
doc nchoosek
nchoosek(3,1)
nchoosek(4,2)
nchoosek(3,2)
edit nchoosek
nchoosk(1:4,2)
nchoosek(1:4,2)
comb2 = nchoosek(1:4,2)
comb2ext = [zeros(6,1), comb2, ones(6,1)*5]
comb2ext(:,2:4) - comb2ext(:,1:3)
comb2ext(:,2:4) - comb2ext(:,1:3) - 1
oc_comb2 = comb2ext(:,2:4) - comb2ext(:,1:3) - 1
kron_index([3 3], indices)
kron_index([3 3], oc_comb2)
kron_index([3 3], oc_comb2+1)
load HTLSs_harmonic
whos
idouble3
kron_index([3 3], oc_comb2+1)
[excitations, indices] = excitation_comb(2, 3)
kron_index([3 3], indices)
excitation_kroni(2, [3 3])
excitation_kroni(1, [3 3])
isingle
isingle3
excitation_kroni(0, [3 3])
[excitations, indices] = excitation_comb(1, 3)
[excitations, indices] = excitation_comb(0, 3)
[excitations, indices] = excitation_comb(2, 3)
doc sparse
a
s3
whos
a3
adag3
Manhar = diag([0 0 -10])
a3s = sparse(a3);
adag3s = sparse(adag3);
Manhars = sparse(diag([0 0 -10]))
whos
H3harc1usIg
H3harcus10nl
H3harcs10nl = multi_kron(eye(3), Hanhar, eye(3)) + multi_kron(eye(3), eye(3), Hanhar) + multi_kron(a3s, adag3s, eye(3)) + multi_kron(adag3s, a3s, eye(3)) + multi_kron(a3s,  eye(3), adag3s) + multi_kron(adag3s,  eye(3), a3s)
H3harcs10nl = multi_kron(eye(3), Manhar, eye(3)) + multi_kron(eye(3), eye(3), Manhar) + multi_kron(a3s, adag3s, eye(3)) + multi_kron(adag3s, a3s, eye(3)) + multi_kron(a3s,  eye(3), adag3s) + multi_kron(adag3s,  eye(3), a3s)
{eye(3), Manhar, eye(3)}
eye(2) + sparse(2,2)
I3s = sparse(eye(3))
H3harcs10nl = multi_kron({I3s, Manhar, I3s}) + multi_kron({I3s, I3s, Manhar}) + multi_kron({a3s, adag3s, I3s}) + multi_kron({adag3s, a3s, I3s}) + multi_kron({a3s, I3s, adag3s}) + multi_kron({adag3s, I3s, a3s})
qubit_excitationsH(H3harcs10nl, [3 3])
H3harcus10nl - ans
%-- 30/11/2020 15:32 --%
load HTLSs_harmonic
Manhar = diag([0 0 -10])
a3s = sparse(a3);
adag3s = sparse(adag3);
Manhars = sparse(diag([0 0 -10]))
I3s = sparse(eye(3))
H3harcs10nl = multi_kron({I3s, Manhar, I3s}) + multi_kron({I3s, I3s, Manhar}) + multi_kron({a3s, adag3s, I3s}) + multi_kron({adag3s, a3s, I3s}) + multi_kron({a3s, I3s, adag3s}) + multi_kron({adag3s, I3s, a3s});
H3LSe_1hus, H3LSe_2hus
H3LSe_1hs
whos
H03e
H3LSe_1hs = multi_kron({I3s, H03e, I3s})
H3LSe_1hus - qubit_excitationsH(H3LSe_1hs, [3 3])
H3LSe_2hs = multi_kron({I3s, I3s, H03e})
H3LSe_2hus - qubit_excitationsH(H3LSe_2hs, [3 3])
save HTLSs_harmonic
H03es = sparse(H03e)
multi_kron({I3s, H03e, I3s})
multi_kron({I3s, I3s, H03e})
generateH2modes
generateH2mode
H2modes
size(H2modes)
size(H2modesu)
H2modesu
isingle2modes = excitation_kroni(1, [3 3 3])
idouble2modes = excitation_kroni(2, [3 3 3])
testTLSmuSG
clear U mniter matvecs
save HTLSs_harmonic
clear all
testTLSmuSG
tic
[time, URK] = ode45(RKfun, [0 T/2 T], ui, options);
toc
toc
testTLSmuSG
load HTLSs_harmonic
clear all
load HTLSs_harmonic
eig(H2modes)
target3LSu
[u0_3 target3LSu]
[excitations2modes1, indices2modes1] = excitation_comb(1, 4)
isingle2modes = excitation_kroni(1, [3 3 3])
[excitations2modes2, indices2modes2] = excitation_comb(2, 4)
idouble2modes = excitation_kroni(2, [3 3 3])
[u0_3 target3LSu]
u02modes = zeros(19,1);
u02modes([1 2 7 11])=1
target2modes(1 11) = 1;
target2modes = zeros(19, 1);
target2modes(1, 11) = 1;
target2modes(3, 6) = 1i
target2modes
target2modes = zeros(19, 1);
target2modes([1, 11]) = 1;
target2modes([3, 6]) = 1i
[u02modes target2modes]
[fieldt2m, fieldw2m, psi2m, relE2m, conv2m, niter2m, mallniterc2m, Jterms2m, maxgrad2m, alpha2m, invHess2m] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-25 105], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
1-Jtermsg4c.Jmax
1-Jmax_overlapg4c
eig(H2modes)
[fieldt2m, fieldw2m, psi2m, relE2m, conv2m, niter2m, mallniterc2m, Jterms2m, maxgrad2m, alpha2m, invHess2m] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-25 105], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 7, 1e-6);
[fieldt2m, fieldw2m, psi2m, relE2m, conv2m, niter2m, mallniterc2m, Jterms2m, maxgrad2m, alpha2m, invHess2m] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-25 105], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
figure
plot(0:0.075:6, fieldtg4c)
hold on
plot(0:0.075:6, fieldt2m)
1-Jtermsg2m.Jmax
1-Jterms2m.Jmax
[fieldt2m1, fieldw2m1, psi2m1, relE2m1, conv2m1, niter2m1, mallniterc2m1, Jterms2m1, maxgrad2m1, alpha2m1, invHess2m1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-25 105], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m1.Jmax
doc divdiff
divdiff
figure
plot(0:0.075:6, fieldt2m1)
eig(H2modes - 6*(H2modes1c+H2modes2c))
eig(H2modes + 6*(H2modes1c+H2modes2c))
[fieldt2m1, fieldw2m1, psi2m1, relE2m1, conv2m1, niter2m1, mallniterc2m1, Jterms2m1, maxgrad2m1, alpha2m1, invHess2m1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-50 120], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m1.Jmax
plot(0:0.075:6, fieldt2m1)
[excitations2modes2, indices2modes2] = excitation_comb(2, 4)
i2harc = [5 9 16:19]
figure
plot(0:0.075:6, conj(psi2m).*psi2m)
plot(0:0.075:6, conj(psi2m(5,:)).*psi2m(5,:))
plot(0:0.075:6, conj(psi2m1(5,:)).*psi2m1(5,:))
plot(0:0.075:6, conj(psi2m1(9,:)).*psi2m1(9,:))
plot(0:0.075:6, conj(psi2m1(16:19,:)).*psi2m1(16:19,:))
[fieldt2m1, fieldw2m1, psi2m1, relE2m1, conv2m1, niter2m1, mallniterc2m1, Jterms2m1, maxgrad2m1, alpha2m1, invHess2m1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-50 120], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m1.Jmax
[fieldt2m1a, fieldw2m1a, psi2m1a, relE2m1a, conv2m1a, niter2m1a, mallniterc2m1a, Jterms2m1a, maxgrad2m1a, alpha2m1a, invHess2m1a] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-50 120], fieldw2m1, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m1a.Jmax
plot(0:0.075:6, fieldt2m1a)
plot(0:0.075:6, conj(psi2m(9,:)).*psi2m(9,:))
plot(0:0.075:6, conj(psi2m(16:19,:)).*psi2m(16:19,:))
options_gate3a = optionsOCqn(1e-4, 0);
save HTLSs_harmonic
options_gate3a.f_max_alpha = get_f_max_alphaOCf_multiE(0.2/0.03, 0.075, 6, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
[~, ~, psi2mg, ~, conv2mg, ~, mallniterc2mg, Jterms2mg, ~, ~, ~] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-50 120], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
[~, ~, psi2mg, ~, conv2mg, ~, mallniterc2mg, Jterms2mg, ~, ~, ~] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes, 2, fcouplingOp2modes, [-50 120], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3a, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mg.Jmax
figure
plot(0:0.075:6, conj(psi2mg(16:19,:)).*psi2mg(16:19,:))
plot(0:0.075:6, conj(psi2m(16:19,:)).*psi2m(16:19,:))
plot(0:0.075:6, conj(psi2m(5,:)).*psi2m(5,:))
plot(0:0.075:6, conj(psi2mg(5,:)).*psi2mg(5,:))
plot(0:0.075:6, conj(psi2m(5,:)).*psi2m(5,:))
plot(0:0.075:6, conj(psi2m(9,:)).*psi2m(9,:))
plot(0:0.075:6, conj(psi2mg(9,:)).*psi2mg(9,:))
plot(0:0.075:6, conj(psi2m(9,:)).*psi2m(9,:))
plot(0:0.075:6, conj(psi2mg(9,:)).*psi2mg(9,:))
plot(0:0.075:6, conj(psi2m(9,:)).*psi2m(9,:))
plot(0:0.075:6, conj(psi2m(5,:)).*psi2m(5,:))
plot(0:0.075:6, conj(psi2mg(5,:)).*psi2mg(5,:))
plot(0:0.075:6, conj(psi2m(16:19,:)).*psi2m(16:19,:))
plot(0:0.075:6, conj(psi2mg(16:19,:)).*psi2mg(16:19,:))
H2modesu
save HTLSs_harmonic
[Jmax_overlap2mg, phases2mg] = Uoverlap_gate(psi2mg(:, end), target2modes, [1 5 9])
[Jmax_overlap2m, phases2m] = Uoverlap_gate(psi2m(:, end), target2modes, [1 5 9])
1-Jmax_overlap4m
1-Jmax_overlap2m
1-Jterms2mg.Jmax
1-Jterms2m.Jmax
1-Jmax_overlap2mg
[6 6; 9 9]
H2modesu([6 6; 9 9])
9
H2modesu([5 5; 9 9])
isingle2ndmode = [5, 9, 16:18];
is_isingle = false(19, 1);
is_isingle(isingle2ndmode) = true;
is_isingle_diag = diag(is_isingle);
is_isingle_diag
H2modesu_i = H2modesu;
omega21 = 20;
H2modesu_i(is_isingle_diag) = omega21;
H2modesu_i(19, 19) = 2*omega21;
H2modesu_i
test_2mode_omega21
ylabel('$\log_{10}(infidelity)$',  'interpreter', 'latex')
xlabel('$\Delta \nu/g$', 'interpreter', 'latex')
hold on
plot(0:0.075:6, fieldt2m)
plot((0:0.075:6)/(2*pi), fieldt2m)
hold on
plot((0:0.075:6)/(2*pi), fieldt2m)
H2modesu10 = H2modesu;
H2modesu10(is_isingle_diag) = 20;
H2modesu10(is_isingle_diag) = 10;
H2modesu10(19,19) = 20;
eig(H2modesu10 + 6*(H2modes1cu+H2modes2cu))
eig(H2modesu10 - 6*(H2modes1cu+H2modes2cu))
Hoperations2mode10 = Hmats2Hops2(H2modesu10, H2modes1cu, H2modes2cu);
H2modesu10
diag(H2modesu10)
[fieldt2m10, fieldw2m10, psi2m10, relE2m10, conv2m10, niter2m10, mallniterc2m10, Jterms2m10, maxgrad2m10, alpha2m10, invHess2m10] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modes10, 2, fcouplingOp2modes, [-25 25], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
[fieldt2m10, fieldw2m10, psi2m10, relE2m10, conv2m10, niter2m10, mallniterc2m10, Jterms2m10, maxgrad2m10, alpha2m10, invHess2m10] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2mode10, 2, fcouplingOp2modes, [-25 25], fieldwg4c, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m10.Jmax
figure
plot(0:0.075:6, fieldt2m)
hold on
plot(0:0.075:6, fieldt2m10)
clf
[fieldt2m10, fieldw2m10, psi2m10, relE2m10, conv2m10, niter2m10, mallniterc2m10, Jterms2m10, maxgrad2m10, alpha2m10, invHess2m10] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2mode10, 2, fcouplingOp2modes, [-25 25], fieldwg4c, @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m10.Jmax
figure
plot(0:0.075:6, fieldt2m10)
[fieldt2m10, fieldw2m10, psi2m10, relE2m10, conv2m10, niter2m10, mallniterc2m10, Jterms2m10, maxgrad2m10, alpha2m10, invHess2m10] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2mode10, 2, fcouplingOp2modes, [-25 25], fieldwg4c, @(w)20/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
plot(0:0.075:6, fieldt2m10)
[fieldt2m10a, fieldw2m10a, psi2m10a, relE2m10a, conv2m10a, niter2m10a, mallniterc2m10a, Jterms2m10a, maxgrad2m10a, alpha2m10a, invHess2m10a] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2mode10, 2, fcouplingOp2modes, [-25 25], fieldw2m10, @(w)20/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
[fieldt2m10a, fieldw2m10a, psi2m10a, relE2m10a, conv2m10a, niter2m10a, mallniterc2m10a, Jterms2m10a, maxgrad2m10a, alpha2m10a, invHess2m10a] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2mode10, 2, fcouplingOp2modes, [-25 25], fieldw2m10, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m10.Jmax
1-Jterms2m10a.Jmax
figure
plot(0:0.075:6, fieldt2m10a)
[fieldt2m10a, fieldw2m10a, psi2m10a, relE2m10a, conv2m10a, niter2m10a, mallniterc2m10a, Jterms2m10a, maxgrad2m10a, alpha2m10a, invHess2m10a] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2mode10, 2, fcouplingOp2modes, [-25 25], fieldw2m10, @(w)50/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 6, 0.075, 9, 9, 1e-6);
1-Jterms2m10a.Jmax
figure
plot(0:0.075:6, fieldt2m10a)
spdiags([0, 0, -50], 1, 1)
doc spdiags
spdiags([0, 0, -50], 0, 3, 3)
spdiags([0, 0, -50], 3, 3, 3)
spdiags([0, 0, -50].', 0, 3, 3)
spdiags([0; 0; -10], 0, 3, 3);
spdiags([0; 0; -10], 0, 3, 3)
a
spdiags(sqrt((1:2).'), 1, 3, 3)
spdiags(sqrt((1:3).'), 1, 3, 3)
spdiags(sqrt((0:2).'), 1, 3, 3)
ans - a3s
spdiags(sqrt((0:2).'), -1, 3, 3)
spdiags(sqrt((1:3).'), -1, 3, 3)
adag3-ans
H03e
spdiags(sqrt((0:2).'), 0, 3, 3)
spdiags((0:2).', 0, 3, 3)
[Hcheck, H1check, H2check, Hucheck, H1ucheck, H2ucheck, Hoperations2modescheck, fcouplingOp2modescheck] = generateH2mode_gen(10, 1, 1, 50);
Hcheck - H2modes
H1check - H2modes1c
H2check - H2modes2c
Hucheck - H2modesu
H1ucheck - H2modes1cu
H2ucheck - H2modes2cu
clear Hcheck H1check H2check Hucheck H1ucheck H2ucheck Hoperations2modescheck fcouplingOp2modescheck
[H2modesa, H1_2modesa, H2_2modesa, Hu2modesa, H1u2modesa, H2u2modesa, Hoperations2modesa, fcouplingOp2modesa] = generateH2mode_gen(0.2/g2GHz, 0.02236/g2GHz, 1, Delta_nuGHz/g2GHz)
g2GHz = 0.03162
Delta_nuGHz = 10.3923 - 9.918
Delta_nuGHz/g2GHz
0.2/g2GHz
0.02236/g2GHz
[H2modesa, H1_2modesa, H2_2modesa, Hu2modesa, H1u2modesa, H2u2modesa, Hoperations2modesa, fcouplingOp2modesa] = generateH2mode_gen(0.2/g2GHz, 0.02236/g2GHz, 1, Delta_nuGHz/g2GHz);
H2modesa
[H2modes H2modesa]
H2modesa
H2modes(29,3)
1/sqrt2
1/sqrt(2)
[H2modesa, H1_2modesa, H2_2modesa, Hu2modesa, H1u2modesa, H2u2modesa, Hoperations2modesa, fcouplingOp2modesa] = generateH2mode_gen(0.2/g2GHz, 1/sqrt(2), 1, 15);
H2modesa
g1GHz*90
g2GHz*90
eig(Hu2modesa + 6*(H1u2modesa+H2u2modesa))
eig(Hu2modesa)
eig(Hu2modesa - 6*(H1u2modesa+H2u2modesa))
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate3, 2.8*2*pi, 0.075, 9, 9, 1e-6);
options_gate4 = optionsOCqn(1e-4, 1e4);
options_gate4.f_max_alpha = get_f_max_alphaOCf_multiE(0.2/0.03, 0.075, 2.8*2*pi, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 0.075, 9, 9, 1e-6);
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 2.8*2*pi/80, 9, 9, 1e-6);
2.8/80
2.8/240
2.8/280
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
options_gate4.f_max_alpha = get_f_max_alphaOCf_multiE(0.2/0.03, 0.01*2*pi, 2.8*2*pi, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)10/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb.Jmax
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb)
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)100/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
max(abs(xmin))
[fieldt2mb, fieldw2mb, psi2mb, relE2mb, conv2mb, niter2mb, mallniterc2mb, Jterms2mb, maxgrad2mb, alpha2mb, invHess2mb] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)50/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb)
plot(0:0.01:2.8, fieldt2mb)
1-Jterms2mb.Jmax
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)100/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb, @(w)100/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb1.Jmax
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)100/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb1.Jmax
figure
plot(0:0.01:2.8, fieldt2mb)
hold on
plot(0:0.01:2.8, fieldt2mb1)
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb1.Jmax
figure
plot(0:0.01:2.8, fieldt2mb1)
mallniterc2mb1
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb, @(w)200/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb1.Jmax
figure
plot(0:0.01:2.8, fieldt2mb1)
[fieldt2mb1, fieldw2mb1, psi2mb1, relE2mb1, conv2mb1, niter2mb1, mallniterc2mb1, Jterms2mb1, maxgrad2mb1, alpha2mb1, invHess2mb1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)100/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb1.Jmax
figure
plot(0:0.01:2.8, fieldt2mb1)
[fieldt2mb2, fieldw2mb2, psi2mb2, relE2mb2, conv2mb2, niter2mb2, mallniterc2mb2, Jterms2mb2, maxgrad2mb2, alpha2mb2, invHess2mb2] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb1, @(w)300/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb2.Jmax
figure
plot(0:0.01:2.8, fieldt2mb2)
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb1, @(w)900/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
figure
plot(0:0.01:2.8, fieldt2mb3)
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb1, @(w)600/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
figure
plot(0:0.01:2.8, fieldt2mb3)
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)600/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)900/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
figure
plot(0:0.01:2.8, fieldt2mb3)
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)650/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)700/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)800/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)850/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)900/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate4, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
options_gate5.f_max_alpha = get_f_max_alphaOCf_multiE(10, 0.01*2*pi, 2.8*2*pi, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
options_gate5 = optionsOCqn(1e-4, 1e4);
options_gate5.f_max_alpha = get_f_max_alphaOCf_multiE(10, 0.01*2*pi, 2.8*2*pi, @(w)320/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), 2);
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb2, @(w)900/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate5, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb3.Jmax
figure
plot(0:0.01:2.8, fieldt2mb3)
[fieldt2mb3, fieldw2mb3, psi2mb3, relE2mb3, conv2mb3, niter2mb3, mallniterc2mb3, Jterms2mb3, maxgrad2mb3, alpha2mb3, invHess2mb3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-0.25/0.03))).*sin(2*pi*w/0.25*0.03), @(w)900/0.03*0.5*(1-tanh(0.6*(w-0.25/0.03))), options_gate5, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
figure
plot(0:0.01:2.8, fieldt2mb3)
figure
plot(0:1/(2*2.8):1/0.02, fieldw2mb)
0.25/0.03
E2modesa = eig(Hu2modesa)
E2modesa = eig(H2modesa)
E2modes_single = eig(H2modesa(2:5,2:5))
E2modes_double = eig(H2modesa(10:19,10:19))
E2modes_single = eig(Hu2modesa(2:5,2:5))
E2modes_double = eig(Hu2modesa(10:19,10:19))
E2modes_double(2:end) - E2modes_double(1:end-1)
[Jmax_overlap2mb, phases2mb] = Uoverlap_gate(psi2mb(:, end), target2modes, [1 5 9])
hold on
plot(0:1/(2*2.8):1/0.02, fieldw2mb1)
plot(0:1/(2*2.8):1/0.02, fieldw2mb2)
plot(0:1/(2*2.8):1/0.02, fieldw2mb3)
plot(0:1/(2*2.8):1/0.02, fieldw2m)
[Jmax_overlap2mb1, phases2mb1] = Uoverlap_gate(psi2mb1(:, end), target2modes, [1 5 9])
[Jmax_overlap2mb2, phases2mb2] = Uoverlap_gate(psi2mb2(:, end), target2modes, [1 5 9])
[Jmax_overlap2mb3, phases2mb3] = Uoverlap_gate(psi2mb3(:, end), target2modes, [1 5 9])
clf
whos
clear invHess2mb invHess2mb1 invHess2mb2 invHess2mb3
save HTLSs_harmonic
plot(0:1/(2*2.8):1/0.02, fieldw2m)
plot(0:1/(2*2.8):1/0.02, fieldw2mb)
H2modesau
Hu2modesa
doc sparse
Hu2modesa+0
[P2ma_single, D2ma_single] = eig(Hu2modesa(2:5,2:5))
[P2ma_single, D2ma_single] = eig(Hu2modesa(2:5,2:5) + 0)
H2modes_singleE =  P2ma_single\Hu2modesa(2:5,2:5)*P2ma_single
H2_2modes_singleE =  P2ma_single\H2u2modesa(2:5,2:5)*P2ma_single
1/0.15
1/0.14
figure
plot(0:1/(2*2.8):1/0.02, fieldw2mb3)
1/0.16
clear H2modes_singleE
E2ma_single =  diag(D2ma_single)
E2ma_single(end) - E2ma_single(1:3)
[P2ma_double, D2ma_double] = eig(Hu2modesa(10:19,10:19) + 0)
H2_2modes_singleE =  P2ma_double\H2u2modesa(10:19,10:19)*P2ma_double
E2ma_double =  diag(D2ma_double)
figure
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3([5, 16:19],:).*psi2mb3([5, 16:19],:)))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3([5, 16:19],:)).*psi2mb3([5, 16:19],:))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3(5,:)).*psi2mb3(5,:))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3(16,:)).*psi2mb3(16,:))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3(17,:)).*psi2mb3(17,:))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3(18,:)).*psi2mb3(18,:))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3(19,:)).*psi2mb3(19,:))
plot(0:1/(2*2.8):1/0.02, conj(psi2mb3([5, 16:19],:)).*psi2mb3([5, 16:19],:))
plot(0:0.01*2*pi:2.8*2*pi, conj(psi2mb3([5, 16:19],:)).*psi2mb3([5, 16:19],:))
plot(0:0.01:2.8, conj(psi2mb3([5, 16:19],:)).*psi2mb3([5, 16:19],:))
figure
plot(0:1/(2*2.8):1/0.02, fieldw2mb)
plot(0:1/(2*2.8):1/0.02, fieldw2mb1)
E2ma_double(2:end) - E2ma_double(2:end-1)
E2ma_double(2:end) - E2ma_double(1:end-1)
E2ma_double(3:end) - E2ma_double(1:end-2)
E2ma_single(2:4) - E2ma_single(1:3)
E2ma_single(3:4) - E2ma_single(1:2)
[fieldt2mb4, fieldw2mb4, psi2mb4, relE2mb4, conv2mb4, niter2mb4, mallniterc2mb4, Jterms2mb4, maxgrad2mb4, alpha2mb4, invHess2mb4] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-6))).*sin(2*pi*w/6), @(w)900/0.03*0.5*(1-tanh(0.6*(w-6))), options_gate5, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
options_gate6 = optionsOCqn(1e-4, 1e4);
options_gate6.f_max_alpha = get_f_max_alphaOCf_multiE(10, 0.01*2*pi, 2.8*2*pi, @(w)0.5*(1-tanh(0.6*(w-6))), 2);
[fieldt2mb4, fieldw2mb4, psi2mb4, relE2mb4, conv2mb4, niter2mb4, mallniterc2mb4, Jterms2mb4, maxgrad2mb4, alpha2mb4, invHess2mb4] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-6))).*sin(2*pi*w/6), @(w)900/0.03*0.5*(1-tanh(0.6*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
figure
1-Jterms2mb4.Jmax
plot(0:1/(2*2.8):1/0.02, fieldw2mb4)
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb4)
[fieldt2mb4, fieldw2mb4, psi2mb4, relE2mb4, conv2mb4, niter2mb4, mallniterc2mb4, Jterms2mb4, maxgrad2mb4, alpha2mb4, invHess2mb4] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(0.6*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(0.6*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb4.Jmax
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb4)
[fieldt2mb5, fieldw2mb5, psi2mb5, relE2mb5, conv2mb5, niter2mb5, mallniterc2mb5, Jterms2mb5, maxgrad2mb5, alpha2mb5, invHess2mb5] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb4, @(w)1e4*0.5*(1-tanh(0.6*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb5.Jmax
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb5)
[fieldt2mb5, fieldw2mb5, psi2mb5, relE2mb5, conv2mb5, niter2mb5, mallniterc2mb5, Jterms2mb5, maxgrad2mb5, alpha2mb5, invHess2mb5] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb4, @(w)1e5*0.5*(1-tanh(0.6*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb5.Jmax
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb5)
1.32/(2*pi)
1/ans
figure
w2ma = 0:1/(2*2.8):1/0.02;
plot(w2ma,
plot(w2ma, 0.5*(1-tanh(0.6*(w2ma-6)))
plot(w2ma, 0.5*(1-tanh(0.6*(w2ma-6))))
plot(w2ma, 0.5*(1-tanh((w2ma-6))))
plot(w2ma, 0.5*(1-tanh(2*(w2ma-6))))
plot(w2ma, 0.5*(1-tanh(2*(w2ma-5))))
figure
plot(w2ma, fieldw2mb5)
options_gate6.f_max_alpha = get_f_max_alphaOCf_multiE(10, 0.01*2*pi, 2.8*2*pi, @(w)0.5*(1-tanh(2*(w-6))), 2);
[fieldt2mb4, fieldw2mb4, psi2mb4, relE2mb4, conv2mb4, niter2mb4, mallniterc2mb4, Jterms2mb4, maxgrad2mb4, alpha2mb4, invHess2mb4] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb4.Jmax
figure
plot(w2ma, fieldw2mb4)
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb4)
[fieldt2mb5, fieldw2mb5, psi2mb5, relE2mb5, conv2mb5, niter2mb5, mallniterc2mb5, Jterms2mb5, maxgrad2mb5, alpha2mb5, invHess2mb5] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
[fieldt2mb5, fieldw2mb5, psi2mb5, relE2mb5, conv2mb5, niter2mb5, mallniterc2mb5, Jterms2mb5, maxgrad2mb5, alpha2mb5, invHess2mb5] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e5*0.5*(1-tanh(2*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-6);
1-Jterms2mb5.Jmax
figure
plot(0:0.01*2*pi:2.8*2*pi, fieldt2mb5)
(10.24-8.796)/(2*pi)
(2*pi)/(10.24-8.796)
figure
plot(w2ma, fieldw2mb5)
figure
plot(w2ma, fieldw2mb)
plot(0:0.01:2.8, fieldt2mb)
plot(0:0.01:2.8, fieldt2mb1)
plot(0:0.01:2.8, fieldt2mb2)
plot(0:0.01:2.8, fieldt2mb3)
figure
plot(0:0.01:2.8, fieldt2mb)
figure
plot(0:0.01:2.8, fieldt2mb1)
1-Jterms2mb1.Jmax
1-Jterms2mb.Jmax
1-Jterms2mb2.Jmax
1-Jterms2mb21.Jmax
1-Jterms2mb1.Jmax
1-Jterms2mb3.Jmax
1-Jterms2mb4.Jmax
1-Jterms2mb5.Jmax
figure
plot(0:0.01:2.8, fieldt2mb5)
xlabel('$g_2t$', 'interpreter', 'latex')
ylabel('$f_n(g_2t)/g_2$',  'interpreter', 'latex')
hold on
plot(0:0.01:2.8, fieldt2mb(2,:))
plot(0:0.01:2.8, fieldt2mb(2,:), 'g')
1-Jterms2mb.Jmax
plot(0:0.01:2.8, fieldt2mb2)
xlabel('$g_2t$', 'interpreter', 'latex')
ylabel('$f_n(g_2t)/g_2$',  'interpreter', 'latex')
doc legend
legend('$n=1$', '$n=2$', 'interpreter', 'latex')
legend({'$n=1$', '$n=2$'}, 'interpreter', 'latex')
1-Jterms2mb.Jmax
1-Jterms2mb2.Jmax
plot(0:0.01:2.8, fieldt2mb5)
xlabel('$g_2t$', 'interpreter', 'latex')
ylabel('$f_n(g_2t)/g_2$',  'interpreter', 'latex')
legend({'$n=1$', '$n=2$'}, 'interpreter', 'latex')
1-Jterms2mb5.Jmax
options_gate7 = optionsOCqn(1e-6, 1e4);
options_gate7.f_max_alpha = get_f_max_alphaOCf_multiE(10, 0.01*2*pi, 2.8*2*pi, @(w)0.5*(1-tanh(2*(w-6))), 2);
options_gate7.invHess0 = invHess2mb5;
[fieldt2mb6, fieldw2mb6, psi2mb6, relE2mb6, conv2mb6, niter2mb6, mallniterc2mb6, Jterms2mb6, maxgrad2mb6, alpha2mb6, invHess2mb6] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modesa, 2, fcouplingOp2modesa, [-25 35], fieldw2mb5, @(w)1e5*0.5*(1-tanh(2*(w-6))), options_gate6, 2.8*2*pi, 1e-2*2*pi, 9, 9, 1e-8);
clear invHess2mb6
save HTLSs_harmonic
whos
figure
plot(t, omega_of_t)
size(t)
t_arbel = [0 t];
omega_of_t = [omega_of_t(end), omega_of_t];
plot(t, omega_of_t)
plot(t_arbel, omega_of_t)
figure
plot(t, t_arbel(2:end)-t_arbel(1:end-1))
plot(t_arbel*g2GHz, omega_of_t)
figure
field_arbel_nu1 = (omega_of_t-10.3923)/g2GHz;
field_arbel_nu2 = (omega_of_t-9.9182)/g2GHz;
figure
plot(t_arbel, field_arbel_nu1)
g2GHz
figure
plot(t_arbel, field_arbel_nu2)
t_arbelg = t_arbel*g2GHz;
plot(t_arbelg, field_arbel_nu1)
plot(t_arbelg, field_arbel_nu2)
xlabel('$g_2t$', 'interpreter', 'latex')
ylabel('$f(g_2t)/g_2$',  'interpreter', 'latex')
figure
plot(t_arbel, field_arbel_nu2)
plot(t_arbel, field_arbel_nu1)
(2.406-0.519)
12/(2*pi)
%-- 08/12/2020 22:43 --%
load HTLSs_harmonic
g1GHzo = 0.0242;
g2GHzo = 0.0291;
[H2modeso, H1_2modeso, H2_2modeso, Hu2modeso, H1u2modeso, H2u2modeso, Hoperations2modeso, fcouplingOp2modeso] = generateH2mode_open(0.2/g1GHzo, 1, g2GHzo/g1GHzo, Delta_nuo);
Delta_nuo = (9.9923 - 10.2876)/g1GHzo
[H2modeso, H1_2modeso, H2_2modeso, Hu2modeso, H1u2modeso, H2u2modeso, Hoperations2modeso, fcouplingOp2modeso] = generateH2mode_open(0.2/g1GHzo, 1, g2GHzo/g1GHzo, Delta_nuo);
options_gate8 = optionsOCqn(1e-4, 1e4);
options_gate8.f_max_alpha = get_f_max_alphaOCf_multiE(6, 0.075, 12, @(w)0.5*(1-tanh(2*(w-6))), 2);
[fieldt2mbo, fieldw2mbo, psi2mbo relE2mbo conv2mbo, niter2mbo, mallniterc2mbo, Jterms2mbo, maxgrad2mbo, alpha2mbo, invHess2mbo] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate8, 12, 0.075, 9, 9, 1e-6);
load HTLSs_harmonic
[fieldt2mbo, fieldw2mbo, psi2mbo relE2mbo conv2mbo, niter2mbo, mallniterc2mbo, Jterms2mbo, maxgrad2mbo, alpha2mbo, invHess2mbo] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate8, 12, 0.075, 9, 9, 1e-6);
1-Jterms2mbo.Jmax
figure
plot(0:0.075:8.55, fieldt2mbo)
plot(0:0.075:12, fieldt2mbo)
[fieldt2mbo1, fieldw2mbo1, psi2mbo1, relE2mbo1, conv2mbo1, niter2mbo1, mallniterc2mbo1, Jterms2mbo1, maxgrad2mbo1, alpha2mbo1, invHess2mbo1] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbo, @(w)1e5*0.5*(1-tanh(2*(w-6))), options_gate8, 12, 0.075, 9, 9, 1e-6);
1-Jterms2mbo1.Jmax
figure
plot(0:0.075:12, fieldt2mbo1)
[fieldt2mbo2, fieldw2mbo2, psi2mbo2, relE2mbo2, conv2mbo2, niter2mbo2, mallniterc2mbo2, Jterms2mbo2, maxgrad2mbo2, alpha2mbo2, invHess2mbo2] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbo, @(w)1e4*0.5*(1-tanh(2*(w-6))), options_gate8, 12, 0.075, 9, 9, 1e-6);
1-Jterms2mbo2.Jmax
figure
plot(0:0.075:12, fieldt2mbo2)
1-Jterms2mbo1.Jmax
whos
options_gate9 = optionsOCqn(1e-4, 1e4);
options_gate9.f_max_alpha = get_f_max_alphaOCf_multiE(6, 0.075, 6, @(w)0.5*(1-tanh(2*(w-6))), 2);
[fieldt2mbo3, fieldw2mbo3, psi2mbo3, relE2mbo3, conv2mbo3, niter2mbo3, mallniterc2mbo3, Jterms2mbo3, maxgrad2mbo3, alpha2mbo3, invHess2mbo3] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e2*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbo3.Jmax
[fieldt2mbo3a, fieldw2mbo3a, psi2mbo3a, relE2mbo3a, conv2mbo3a, niter2mbo3a, mallniterc2mbo3a, Jterms2mbo3a, maxgrad2mbo3a, alpha2mbo3a, invHess2mbo3a] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbo3, @(w)1e4*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbo3a.Jmax
figure
plot(0:0.075:6, fieldt2mbo3a)
[fieldt2mbo3b, fieldw2mbo3b, psi2mbo3b, relE2mbo3b, conv2mbo3b, niter2mbo3b, mallniterc2mbo3b, Jterms2mbo3b, maxgrad2mbo3b, alpha2mbo3b, invHess2mbo3b] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbo3, @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbo3b.Jmax
figure
plot(0:0.075:6, fieldt2mbo3b)
[fieldt2mbo3c, fieldw2mbo3c, psi2mbo3c, relE2mbo3c, conv2mbo3c, niter2mbo3c, mallniterc2mbo3c, Jterms2mbo3c, maxgrad2mbo3c, alpha2mbo3c, invHess2mbo3c] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbo3, @(w)5e2*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbo3b.Jmax
figure
plot(0:0.075:6, fieldt2mbo3c)
1-Jterms2mbo3c.Jmax
save HTLSs_harmonic
target2mswap = zeros(19,1);
target2mswap([1, 3 6 11]) = 1;
[target2modes target2mswap]
[fieldt2mbos, fieldw2mbos, psi2mbos, relE2mbos, conv2mbos, niter2mbos, mallniterc2mbos, Jterms2mbos, maxgrad2mbos, alpha2mbos, invHess2mbos] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate8, 12, 0.075, 9, 9, 1e-6);
1-Jterms2mbos.Jmax
figure
plot(0:0.075:12, fieldt2mbos)
[fieldt2mbos1, fieldw2mbos1, psi2mbos1, relE2mbos1, conv2mbos1, niter2mbos1, mallniterc2mbos1, Jterms2mbos1, maxgrad2mbos1, alpha2mbos1, invHess2mbos1] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos, @(w)1e4*0.5*(1-tanh(2*(w-6))), options_gate8, 12, 0.075, 9, 9, 1e-6);
1-Jterms2mbos1.Jmax
figure
plot(0:0.075:12, fieldt2mbos1)
options_gate10 = optionsOCqn(1e-4, 1e4);
options_gate10.f_max_alpha = get_f_max_alphaOCf_multiE(6, 0.075, 8.55, @(w)0.5*(1-tanh(2*(w-6))), 2);
[fieldt2mbo4, fieldw2mbo4, psi2mbo4, relE2mbo4, conv2mbo4, niter2mbo4, mallniterc2mbo4, Jterms2mbo4, maxgrad2mbo4, alpha2mbo4, invHess2mbo4] = OClimf_gate(u02modes, target2modes, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate10, 8.55, 0.075, 9, 9, 1e-6);
1-Jterms2mbo4.Jmax
figure
plot(0:0.075:8.55, fieldt2mbo4)
hold on
plot(0:0.075:8.55, fieldt2mbo4(:, end:-1:1))
max(abs(fieldt2mbo4-fieldt2mbo4(:, end:-1:1)), [], 2)
max(abs(fieldt2mbo4-fieldt2mbo4([2 1], end:-1:1)), [], 2)
[Jmax_overlap2mbo4, phases2mbo4] = Uoverlap_gate(psi2mbo4(:, end), target2modes, [1 5 9])
figure
plot(0:0.075:8.55, conj(psi2mbo4).*psi2mbo4)
plot(0:0.075:8.55, conj(psi2mbo4(5,:)).*psi2mbo4(5,:))
plot(0:0.075:8.55, conj(psi2mbo4(16:19,:)).*psi2mbo4(16:19,:))
save HTLSs_harmonic
1-Jterms2mbo.Jmax
1-Jterms2mbo1.Jmax
1-Jterms2mbo2.Jmax
1-Jterms2mbos.Jmax
1-Jterms2mbos1.Jmax
1-Jterms2mbo4.Jmax
1-Jterms2mbo3c.Jmax
1-Jterms2mbo3.Jmax
1-Jterms2mbo3a.Jmax
1-Jterms2mbo3b.Jmax
figure
plot(0:0.075:6, fieldt2mbo3a)
plot(0:0.075:6, fieldt2mbo3b)
plot(0:0.075:6, fieldt2mbo3a)
1-Jterms2mbos.Jmax
1-Jterms2mbos1.Jmax
figure
plot(0:0.075:6, fieldt2mbos1)
plot(0:0.075:12, fieldt2mbos1)
plot(0:0.075:12, fieldt2mbos2)
plot(0:0.075:12, fieldt2mbos1a)
[fieldt2mbos2, fieldw2mbos2, psi2mbos2, relE2mbos2, conv2mbos2, niter2mbos2, mallniterc2mbos2, Jterms2mbos2, maxgrad2mbos2, alpha2mbos2, invHess2mbos2] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
figure
plot(0:0.075:6, fieldt2mbos2)
1-Jterms2mbos2.Jmax
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)1e4*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbos3.Jmax
figure
plot(0:0.075:6, fieldt2mbos3)
[fieldt2mbos2, fieldw2mbos2, psi2mbos2, relE2mbos2, conv2mbos2, niter2mbos2, mallniterc2mbos2, Jterms2mbos2, maxgrad2mbos2, alpha2mbos2, invHess2mbos2] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], @(w)0.05*0.5*(1-tanh(2*(w-6))).*sin(2*pi*w/6), @(w)1e2*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbos2.Jmax
figure
plot(0:0.075:6, fieldt2mbos2)
max(abs(fieldt2mbos2-fieldt2mbos2([2 1], end:-1:1)), [], 2)
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)1e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
figure
plot(0:0.075:6, fieldt2mbos3)
1-Jterms2mbos3.Jmax
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)1e4*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbos3.Jmax
figure
plot(0:0.075:6, fieldt2mbos3)
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)5e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbos3.Jmax
figure
plot(0:0.075:6, fieldt2mbos3)
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)2e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbos3.Jmax
figure
plot(0:0.075:6, fieldt2mbos3)
[fieldt2mbos3, fieldw2mbos3, psi2mbos3, relE2mbos3, conv2mbos3, niter2mbos3, mallniterc2mbos3, Jterms2mbos3, maxgrad2mbos3, alpha2mbos3, invHess2mbos3] = OClimf_gate(u02modes, target2mswap, [1 5 9], Hoperations2modeso, 2, fcouplingOp2modeso, [-30 30], fieldw2mbos2, @(w)2.5e3*0.5*(1-tanh(2*(w-6))), options_gate9, 6, 0.075, 9, 9, 1e-6);
1-Jterms2mbos3.Jmax
plot(0:0.075:6, fieldt2mbos3)
save HTLSs_harmonic