[P, D] = diag([0 1i; 1i 0])
[P, D] = eig([0 1i; 1i 0])
sqrtD= sqrt(D)
P*sqrtD
sqrtD.^2
P*sqrtD\P
ans.^2
P*sqrtD*inv(P)
P*sqrtD*P'
ans.^2
sqrtiswap = P*sqrtD*P'
sqrtiswap^2
speye(4)
omegac0 = 7445/71.4;
omega1 = 5889.9/71.4;
omega2 = 5031.1/71.4;
alpha1 = 324/71.4;
alpha2 = 235/71.4;
g1 = 100/71.4;
g2 = 1;
n = (0:3).';
N_HO = spdiags(n, 0, 4, 4);
Manhar = spdiags([0; 0; -1; -3], 0, 4, 4);
Mcoup = spdiags(sqrt([(n + 1).', n]), [-1, 1], 4, 4);
I4 = speye(4);
Mcoup = spdiags(sqrt([(n + 1).', n]), [-1, 1], 4, 4);
I4 = speye(4);
Mcoup = spdiags(sqrt([(n + 1), n]), [-1, 1], 4, 4);
I4 = speye(4);
Mcoup
multi_kron({omegac*N_HO + alphac*Manhar, I4, I4})
alphac = 230/71.4;
clear omegac
omegac0 = 7445/71.4;
multi_kron({omegac*N_HO + alphac*Manhar, I4, I4})
multi_kron({omegac0*N_HO + alphac*Manhar, I4, I4})
H0 = multi_kron({omegac0*N_HO + alphac*Manhar, I4, I4}) + multi_kron({I4, omega1*N_HO + alpha1*Manhar, I4})...
+ multi_kron({I4, I4, omega2*N_HO + alpha2*Manhar})...
+ g1*multi_kron({Mcoup, Mcoup, I4}) + g2*multi_kron({Mcoup, I4, Mcoup});
eig(H0)
doc delta
T = 0.1142*71.4*2*pi
sigma = 8.3e-3*71.4*2*pi
T/(2*pi)
sigma/(2*pi)
T = 0.1142*100*2*pi
sigma = 8.3e-3*100*2*pi
sigma/(2*pi)
busHamiltonian2q
eig(H0)
t=0:T/1e3:T;
plot(t, 0.5*(1-tanh(t-3*sigma)).*(t-3*sigma) + 0.5*(1 + tanh(t-T+3*sigma)).*(t-T+3*sigma))
3*sigma
plot(t, 0.5*(1-tanh(t-3.2*sigma)).*(t-3*sigma) + 0.5*(1 + tanh(t-T+3.2*sigma)).*(t-T+3*sigma))
plot(t, 0.5*(1-tanh(t-2.8*sigma)).*(t-3*sigma) + 0.5*(1 + tanh(t-T+2.8*sigma)).*(t-T+3*sigma))
plot(t, 0.5*(x-log(cosh(t-3*sigma))) + 0.5*(x + log(cosh(t-T+3*sigma))))
plot(t, 0.5*(t-log(cosh(t-3*sigma))) + 0.5*(t + log(cosh(t-T+3*sigma))))
0.5*(T/2-log(cosh(T/2-3*sigma))) + 0.5*(T/2 + log(cosh(T/2-T+3*sigma)))
plot(t, t + 0.5*log(cosh(t-T+3*sigma)/cosh(t-3*sigma)))
plot(t, t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma)))
T/2 + 0.5*log(cosh(T/2-T+3*sigma)./cosh(T/2-3*sigma))
midf = T/2 + 0.5*log(cosh(T/2-T+3*sigma)./cosh(T/2-3*sigma))
plot(t, t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma))-midf)
3*sigma
2*pi/omegac0
omegac0
T
tg = t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma))-midf;
figure
plot(t, exp(-tg.^2/(2*sigma^2)))
envc = exp(-tg.^2/(2*sigma^2));
envdc = zeros(1, 1001);
T/(3*sigma)
(3*sigma)/T*1e3
218*T/1e3
3*sigma
envdc(1:219) = exp(-t.^2/(2*sigma^2));
envdc(1:219) = exp(-t(1:219).^2/(2*sigma^2));
1001-218
envdc(783:1001) = exp(-(t(783:1001)-T+3*sigma).^2/(2*sigma^2));
envdc(1:219) = exp(-(t(1:219)-3*sigma).^2/(2*sigma^2));
envdc(220:782) = 1;
hold on
plot(t, envdc)
figure
plot(t, envdc - envc)
midf -T/2
midf
midf - T/2
envelope = @(t) exp(-(t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma))-T/2).^2/(2*sigma^2))
max(abs(envc-envelope(t)))
gs4 = sparse(4, 1);
gs4(1) = 1;
fundamental4 = sparse(4, 1);
fundamental4(2) = 1;
down_down = multi_kron({gs4, gs4, gs4});
down_up = multi_kron({gs4, gs4, fundamental4});
up_down = multi_kron({gs4, fundamental4, gs4});
up_up = multi_kron({gs4, fundamental4, fundamental4});
U2q0 = [down_down, down_up, up_down, up_up];
Usq_iswap = [down_down, (down_up + 1i*up_down)/sqrt(2), (1i*down_up + up_down)/sqrt(2), up_up];
ui = sparse(256, 1);
ui([1, 66, 133, 262]) = 1;
max(abs(U2q0(:)-ui))
size(U2q0)
size(U2q0(:))
size(ui)
ui = sparse(256, 1);
ui([1, 66, 133, 198]) = 1;
max(abs(U2q0(:)-ui))
clear Usq_iswap;
Usqrt_iswap = [down_down, (down_up + 1i*up_down)/sqrt(2), (1i*down_up + up_down)/sqrt(2), up_up];
Usqrt_iswap
size(Usqrt_iswap)
U2q0
eig(H0)
eig(H0 - Ncoupler)
T
7*530
tgrid = 0:T/1e3:T;
options = SGdefault_op;
data = SGdata(options);
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 4e3, 9, 7, 1e-5);
busHamiltonian2q
mniter
T/4e3
ans*530
T*530
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 4e3, 9, 7, 1e-5);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
max(abs(ev_domain))
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 7, 1e-5);
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
mniter
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 9, 1e-5);
max(abs(eig(H0u)))
max(abs(control(t)))
max(abs(control(tgrid)))
busHamiltonian2q
mniter
est_errors
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 4e3, 9, 9, 1e-5);
mniter
est_errors
sqnorm(U(:,end))
sqnorm(U(:,end).*Usqrt_iswap(:))
U(1:8,end).*conj(U(1:8,end))
U(:,end).*conj(U(:,end))
[U(:,end).*conj(U(:,end)), (1:256).']
nnz(U(:,end))
nnz(U(:,end).*conj(U(:,end)))
nnz(U2q0(:))
nnz(U(1:64,end).*conj(U(1:64,end)))
nnz(U(65:128,end).*conj(U(65:128,end)))
nnz(U(129:192,end).*conj(U(129:192,end)))
nnz(U(193:256,end).*conj(U(193:256,end)))
reshape(U(:, end).*conj(U(:, end)), [64,4])
[reshape(U(:, end).*conj(U(:, end)), [64,4]), (1:64).']
sqnorm(conj(U(:,end)).*Usqrt_iswap(:))
norm(conj(U(:,end)).*Usqrt_iswap(:))
norm(conj(U(:,end)).*Usqrt_iswap(:))^2
norm(conj(U(65:128,end)).*Usqrt_iswap(:,2))^2
sqnorm(Usqrt_iswap)
norm(conj(U(129:192,end)).*Usqrt_iswap(:,3))^2
abs(Usqrt_iswap'*U)^2
abs(Usqrt_iswap'*U).^2
abs(Usqrt_iswap(:)'*U).^2
sqnorm(Usqrt_iswap(:).*U(:,end))
sqnorm(Usqrt_iswap(:).*U)
figure
plot(t, control(t))
omega_phi
T
U(65:128, end)
busHamiltonian2q
norm(conj(U(129:192,end)).*Usqrt_iswap(:,3))^2
norm(conj(U(:,end)).*Usqrt_iswap(:))^2
figure
plot(t, theta + delta*envelope(t).*cos(omega_phi*t))
N_HO
Manhar
Mcoup
figure
plot(t, U(65:128,:).*conj(U(65:128,:)))
plot(t, U(129:192,:).*conj(U(129:192,:)))
norm(conj(abs(U(:,end))).*abs(Usqrt_iswap(:)))^2
norm(abs(U(:,end).*Usqrt_iswap(:)))^2
sum(abs(U(:,end).*Usqrt_iswap(:)).^2)
abs(Usqrt_iswap(:))
abs(U(:,end))
abs(reshape(U(:,end)[64,4]))
abs(reshape(U(:,end),[64,4]))
abs(reshape(U(:,end),[64,4]).*(Usqrt_iswap))
sum(abs(reshape(U(:,end),[64,4]).*(Usqrt_iswap)))
sum(abs(U(:,end).*Usqrt_iswap(:,end)))
sum(abs(U(:,end).*Usqrt_iswap(:)))
sum(abs(U(:,end).*Usqrt_iswap(:)))/4
busHamiltonian2q
sum(abs(U(:,end).*Usqrt_iswap(:)))/4
figure
plot(tgrid, U(129:192,:).*conj(U(129:192,:)))
62.05/(100*2*pi)
busHamiltonian2q
est_errors
mniter
matvecs
fidelity = sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16;
infidelity = 1 - fidelity;
fidelity
max(history.niter)
40^9/factorial(9)
busHamiltonian2q
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
(0.1/ans)^(1/18)
[Uex, mniter_ex, matvecs_ex, est_errors_ex, history_ex] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 2e4, 9, 13, eps);
est_errors_ex
mniter_ex
[Uex, mniter_ex, matvecs_ex, est_errors_ex, history_ex] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 3e4, 9, 13, eps);
est_errors_ex
excitation_kroni(0, [4, 4])
excitation_kroni(2, [4, 4])
excitation_kroni(4, [4, 4])
excitation_kroni(6, [4, 4])
size(excitation_kroni(6, [4, 4]))
parent_indices = inv_kron_index([4, 4], (1:64))
parent_indices = inv_kron_index([4, 4], (1:64).')
sum(parent_indices, 2)
sum(parent_indices-1, 2)
index2exitations = sum(parent_indices-1, 2);
[reshape(U(:, end).*conj(U(:, end)), [64,4]), index2exitations]
norm(U(:,end)-Uex(:,end))
norm(U(:,end)-Uex(:,end))/norm(U(:,end))
est_errors
tic
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 5, 1e-5);
toc
S2q2full = [down_down, down_up, up_down, up_up];
Usqrt_iswap2q = sparse([1 0])
Usqrt_iswap2q = sparse([1 0 0 0; 0 1/sqrt(2) 1i/sqrt(2) 0; 0 1i/sqrt(2) 1/sqrt(2), 0; 0 0 0 1])
Usqrt_iswap - S2q2full*Usqrt_iswap2q*S2q2full.'
S2q2full*Usqrt_iswap2q*S2q2full.'
size(ans)
size(S2q2full)
size(S2q2full*Usqrt_iswap2q*S2q2full')
size(S2q2full*(Usqrt_iswap2q*S2q2full'))
size((Usqrt_iswap2q*S2q2full'))
size(S2q2full(:, 1:4)*(Usqrt_iswap2q*S2q2full'))
Usqrt_iswap - S2q2full*Usqrt_iswap2q
norm(abs(U(:,end).*Usqrt_iswap(:)))^2
sum(abs(U(:,end).*Usqrt_iswap(:)))/4
busHamiltonian2q
sum(abs(U(:,end).*Usqrt_iswap(:)))/4
figure
plot(tgrid, U(129:192,:).*conj(U(129:192,:)))
sqnorm(Usqrt_iswap(:).*U)
sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16
fidelity
sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16 - fidelity
theta
theta/pi
figure
plot(t, theta + delta*envelope(t).*cos(omega_phi*t))
hold on
plot(t, envelope(t).*(theta + delta*cos(omega_phi*t)))
figure
plot(t, omegac0*(1 - sqrt(cos(pi*(theta + delta*envelope(t).*cos(omega_phi*t))))))
hold on
plot(t, omegac0*(1 - sqrt(cos(pi*envelope(t).*(theta + delta*cos(omega_phi*t))))))
busHamiltonian2q
sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16;
sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16
busHamiltonian2q
sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16
figure
plot(tgrid, U(129:192,:).*conj(U(129:192,:)))
sum(Usqrt_iswap.*reshape(U(:, end), [64 4]).^2
sum(Usqrt_iswap.*reshape(U(:, end), [64 4])).^2
abs(sum(Usqrt_iswap.*reshape(U(:, end), [64 4]))).^2
(sum(abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4]))).^2
sum((sum(abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4])))).^2/16
sum((sum(abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4]))).^2)/4
busHamiltonian2q
sum((sum(abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4])))).^2/16
(sum(abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4]))).^2
plot(tgrid, U(129:192,:).*conj(U(129:192,:)))
plot(tgrid, U(65:128,:).*conj(U(65:128,:)))
reshape(abs(U(:, end)), [64 4]))
reshape(abs(U(:, end)), [64 4])
plot(tgrid, U(65:128,:).*conj(U(65:128,:)))
plot(tgrid, U(129:192,:).*conj(U(129:192,:)))
abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4]))
abs(Usqrt_iswap).*reshape(abs(U(:, end)), [64 4])
abs(Usqrt_iswap)
bus
busHamiltonian2q
fidelity = sum(abs(U(:,end).*Usqrt_iswap(:)))^2/16;
fidelity
U([])
U([1:0])
prod(U([1:0]))
A
A = [1, 2, 3; 4 5 6; 7 8 9]
A([2,3; 4,5])=0
A = [1, 2, 3; 4 5 6; 7 8 9]
A([2,3])=0
A = [1, 2, 3; 4 5 6; 7 8 9]
A([2,4], [3, 5]) = 0
A = [1, 2, 3; 4 5 6; 7 8 9]
A([2,3], [3, 3]) = 0
U2q0test = generateU2q0full([4,4,4]);
U2q0test
A = [1, 2, 3; 4 5 6; 7 8 9]
A([2,3], [2, 3]) = 0
U2q0test = generateU2q0full([4,4,4])
U2q0test-U2q0
U2q0test = generateU2q0full([4,4,3])
size(U2q0test)
clear U2q0test
Usqrt_iswap_2q = sparse([1, 0,          0,          0;
0, 1/sqrt(2),  1i/sqrt(2), 0;
0, 1i/sqrt(2), 1/sqrt(2),  0;
0, 0,          0,          0]);
Usqrt_iswap - U2q0*Usqrt_iswap_2q
U2q0
Usqrt_iswap_2q = sparse([1, 0,          0,          0;
0, 1/sqrt(2),  1i/sqrt(2), 0;
0, 1i/sqrt(2), 1/sqrt(2),  0;
0, 0,          0,          1]);
Usqrt_iswap - U2q0*Usqrt_iswap_2q
Nex = Nexcitations_kron([4, 4], (1:64).');
[reshape(U(:, end).*conj(U(:, end)), [64,4]), Nex]
[reshape(U(:, end).*conj(U(:, end)), [64,4]), index2exitations]
[reshape(U(:, end).*conj(U(:, end)), [64,4]), Nex]
Nex = Nexcitations_kron([4, 4], (1:64).');
[reshape(U(:, end).*conj(U(:, end)), [64,4]), Nex]
index2exitations - Nex
Nex_is_even = (mod(Nex, 2) == 0)
U4end = reshape(U(:,end), [64,4]);
sqnorm(U4end)
sqnorm(U4end(Nex_is_even, 1))
sqnorm(U4end(Nex_is_even, 2))
sqnorm(U4end(Nex_is_even, 3))
sqnorm(U4end(Nex_is_even, 4))
sqnorm(U4end(~Nex_is_even, 2))
sqnorm(U4end(~Nex_is_even, 3))
H0even = H0(Nex_is_even, Nex_is_even)
size(H0even)
max(max(abs(H0even-H0even')))
Nex_is_odd = ~Nex_is_even;
H0odd = H0(Nex_is_odd, Nex_is_odd);
H0odd
H0u_par = sparse(128)
H0u_par = sparse(128, 128)
B = blkdiag(A, rand(2), A)
U2q0test = generateU2q0full([4,4,3])
U2q0test = generateU2q0full([4,4,4])
U2q0test-U2q0
U2q0test = generateU2q0full([4,4,4])
U2q0test-U2q0
H0u_par = sparse(blkdiag(H0even, H0odd, H0odd, H0even));
Ncoupler_even = Ncoupler(Nex_is_even, Nex_is_even);
Ncoupler_even
Ncoupler_odd = Ncoupler(Nex_is_odd, Nex_is_odd);
Ncoupler_odd
Ncoupler_even - Ncoupler_odd
Ncoupler_u_par = kron(I4, Ncoupler_even);
Gop_par = @(u, t, v) -1i*(H0u_par*v - control(t)*Ncoupler_u_par*v);
Gdiff_op_par = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler_u_par*u1);
U2q0par = [U2q0(Nex_is_even),U2q0(Nex_is_odd), U2q0(Nex_is_odd), U2q0(Nex_is_even)]
U2q0par = [U2q0(Nex_is_even, 1),U2q0(Nex_is_odd, 2), U2q0(Nex_is_odd, 3), U2q0(Nex_is_even, 1)]
U2q0par = [U2q0(Nex_is_even, 1),U2q0(Nex_is_odd, 2), U2q0(Nex_is_odd, 3), U2q0(Nex_is_even, 4)]
[Up, mniterp, matvecsp, est_errorsp, historyp] = SemiGlobal1(Gop_par, Gdiff_op_par, 0, [], [-528*1i, 1i], U2q0par(:), tgrid, 8e3, 9, 5, 1e-5);
tic, [Up, mniterp, matvecsp, est_errorsp, historyp] = SemiGlobal1(Gop_par, Gdiff_op_par, 0, [], [-528*1i, 1i], U2q0par(:), tgrid, 8e3, 9, 5, 1e-5); toc
5/13
restore_vec = [Nex_is_even; Nex_is_odd; Nex_is_odd; Nex_is_even];
size(restore_vec)
sqnorm(Up(:,end))
Up_restored(restor_vec, :) = Up;
Up_restored(restore_vec, :) = Up;
size(Up_restore)
size(Up_restored)
size(restore_vec)
size(Up)
Up_restored = zeros(256, 1001);
Up_restored(restore_vec, :) = Up;
size(restore_vec)
size(Up_restored)
max(max(abs(U-Up_restored)))
mniterp
mniter
matvecsp
matvecs
Gop_psi = @(u, t, v) -1i*(H0*v - control(t)*Ncoupler*v);
Gdiff_op_psi = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler*u1);
tic, [psi2, mniter2, matvecs2, est_errors2, history2] = SemiGlobal1(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
3.778*4
tic, [psi2, mniter2, matvecs2, est_errors2, history2] = SemiGlobal1(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
3.479*4
tic, [psi2, mniter2, matvecs2, est_errors2, history2] = SemiGlobal1(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
tic
[U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 5, 1e-5);
tic, [U, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 5, 1e-5); toc
mniter2
tic, for j=1:4, [psij, mniterj, matvecsj, est_errorsj, historyj] = SemiGlobal1(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0(:, j), tgrid, 8e3, 9, 5, 1e-5);end, toc
mniterj
est_errorsj
max(max(abs(U(65:128,:)-psi2)))
max(max(abs(U(193:256,:)-psij)))
clear psij mniterj matvecsj est_errorsj historyj
busHamiltonian2q
max(max(abs(Upar-Up)))
Ucopy = U;
busHamiltonian2q
max(max(abs(U-Ucopy)))
13.5/15
tic, [psip2, mniterp2, matvecsp2, est_errorsp2, historyp2] = SemiGlobal1(Gop_p_psi, Gdiff_op_p_psi, 0, [], [-528*1i, 1i], U2q0par(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
Gop_even = @(u, t, v) -1i*(H0even*v - control(t)*Ncoupler_u*v);
Gdiff_op_even = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler_even*u1);
Gop_odd = @(u, t, v) -1i*(H0odd*v - control(t)*Ncoupler_even*v);
Gop_even = @(u, t, v) -1i*(H0even*v - control(t)*Ncoupler_even*v);
tic, [psip2, mniterp2, matvecsp2, est_errorsp2, historyp2] = SemiGlobal1(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
size(psip2)
max(max(abs(Upar(65:128,:)-psip2)))
max(max(abs(Upar(33:64,:)-psip2)))
tic, [psip4, mniterp4, matvecsp4, est_errorsp4, historyp4] = SemiGlobal1(Gop_even, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 4), tgrid, 8e3, 9, 5, 1e-5);toc
max(max(abs(Upar(97:128,:)-psip4)))
tic, zeros(256, 64e3); toc
tic, B=zeros(256, 64e3); toc
clear B
[psip4, mniterp4, matvecsp4, est_errorsp4, historyp4] = SemiGlobal1(Gop_even, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 4), tgrid, 8e3, 9, 5, 1e-5);
tic, [psi2, mniter2, matvecs2, est_errors2, history2] = SemiGlobal1(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
[psi2, mniter2, matvecs2, est_errors2, history2] = SemiGlobal1(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0(:, 2), tgrid, 8e3, 9, 5, 1e-5);
H0u_full = kron(I4, H0);
Ncoupler_u_full = kron(I4, Ncoupler);
[Upar, mniter, matvecs, est_errors, history] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0par(:), tgrid, 8e3, 9, 5, 1e-5, options, data);
Gop_full = @(u, t, v) -1i*(H0u_full*v - control(t)*Ncoupler_u_full*v);
Gdiff_op_full = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler_u_full*u1);
tic, [Uf, mniterf, matvecsf, est_errorsf, historyf] = SemiGlobal1(Gop_full, Gdiff_op_full, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 5, 1e-5, options, data);toc
[Uf, mniterf, matvecsf, est_errorsf, historyf] = SemiGlobal1(Gop_full, Gdiff_op_full, 0, [], [-528*1i, 1i], U2q0(:), tgrid, 8e3, 9, 5, 1e-5, options, data);
max(abs(control(tgrid)))
omegac0
omegac0-omega1
a4 = spdiags(sqrt(n), 1, 4, 4)
adag4 = spdiags(sqrt(n + 1), -1, 4, 4)
busHamiltonian2qI
size(b_a1)
b_a1(Nex_is_even)
b_a1(Nex_is_odd)
b_a1
b_a1(Nex_is_even, Nx_is_even)
b_a1(Nex_is_even, Nex_is_even)
b_a1(Nex_is_odd, Nex_is_odd)
b_a1(Nex_is_odd, Nex_is_odd)- b_a1(Nex_is_even, Nex_is_even)
b_a2(Nex_is_odd, Nex_is_odd)- b_a2(Nex_is_even, Nex_is_even)
b_a1
b_a2
b_a1(Nex_is_even, Nex_is_even)
b_a2(Nex_is_even, Nex_is_even)
eig(H0I)
max(abs(control(tgrid)))
controlI = @(t) Dc0*(1 - sqrt(cos(pi*(theta + delta*envelope(t).*cos(omega_phi*t)))));
min(abs(control(tgrid)))
Dc0
eig(H0I - 3*Ncoupler_u)
eig(H0I - 3*Ncoupler)
eig(H0I + bdag_a1dag)
eig(full(H0I + bdag_a1dag))
max(ans)
min(eig(full(H0I + bdag_a1dag)))
busHamiltonian2qI
sqnorm(UIpar(:,end))
infidelityI
fidelityI
mniter
est_errorsI
matvecsI
matvecs
mniterI
UI(:,end).*conj(UI(:,end))
reshape(UI(:,end).*conj(UI(:,end)), [64,4])
busHamiltonian2qI
eig(full(H0I + bdag_a1dag))
eig(H0I)
eig(H0I - 3*Ncoupler)
eig(full(H0I + bdag_a1dag))
max(ans)
min(ans)
min(eig(full(H0I + bdag_a1dag)))
busHamiltonian2qI
reshape(UI(:,end).*conj(UI(:,end)), [64,4])
fidelityI
mniterI
Dc0
max(abs(controlI(tgrid)))
D2
omega1
sqnorm(UIpar(:,end))
b_a1even
b_a1odd
b_a2even
b_a2odd
busHamiltonian2qI
fidelityI
busHamiltonian2qI
max(abs(control(tgrid)))
eig(H0I - 13*Ncoupler)
busHamiltonian2qI
fidelityI
busHamiltonian2qI
mniterI
fidelityI
fidelity
reshape(UI(:,end).*conj(UI(:,end)), [64,4])
busHamiltonian2qI
fidelity
fidelityI
fidelity-fidelityI
max(abs(U(:,end) - fMdiag(H0, @(H) -1i*H, T, 1, UI(:,end)))
max(abs(U(:,end) - fMdiag(H0, @(H) -1i*H, T, 1, UI(:,end))))
UIback = fMdiag(H0, @(H) -1i*H, UI(:,end), T, 1);
UIback = fMdiag(H0, @(H) exp(-1i*H), UI(:,end), T, 1);
size(H0)
[P,D] = eig(full(H0));
UIback = P*exp(-1i*D*T)*P\reshape(UI(:,end), [64,4]);
max(abs(U(:,end) - UIback(:))
max(abs(U(:,end) - UIback(:)))
UIback
sqnorm(UIback)
UIback = P*exp(-1i*D*T)*inv(P)*reshape(UI(:,end), [64,4]);
sqnorm(UIback)
diag(D)
UIback = P*exp(-1i*D*T)*P'*reshape(UI(:,end), [64,4]);
sqnorm(UIback)
sqnorm(reshape(UI(:,end), [64,4]))
sqnorm(P)
sqnorm(exp(-1i*D*T))
D
(exp(-1i*D*T))
evH0 = diag(D);
UIback = P*evH0.*P'*reshape(UI(:,end), [64,4]);
sqnorm(reshape(UI(:,end), [64,4]))
sqnorm(UIback)
UIback = P*exp(-1i*evH0*T).*P'*reshape(UI(:,end), [64,4]);
sqnorm(UIback)
UIback = P*(exp(-1i*evH0*T).*(P'*reshape(UI(:,end), [64,4])));
sqnorm(UIback)
max(abs(U(:,end) - UIback(:)))
UIback
abs(UIback)
[PI, DI] = eig(H1har);
H1har = multi_kron({I4, omega1*N_HO, I4});
[PI, DI] = eig((full(H1har));
[PI, DI] = eig((full(H1har)));
evH1har = diag(DI);
UIback = PI*(exp(-1i*evH1har*T).*(PI'*reshape(UI(:,end), [64,4])));
sqnorm(UIback)
max(abs(U(:,end) - UIback(:)))
evH1har
max(abs(U(:,end)) - abs(UIback(:))))
max(abs(U(:,end)) - abs(UIback(:)))
diag(H1har)
UIback = exp(-1i*diag(H1har)*T).*reshape(UI(:,end), [64,4]);
sqnorm(UIback)
max(abs(U(:,end)) - abs(UIback(:)))
max(abs(U(:,end) - UIback(:)))
UIback = exp(1i*diag(H1har)*T).*reshape(UI(:,end), [64,4]);
max(abs(U(:,end) - UIback(:)))
Hint = multi_kron({I4, omega1*N_HO, I4}) + multi_kron({omega1*N_HO, I4, I4}) + multi_kron({I4, I4, omega1*N_HO});
Hint
UIback = exp(1i*diag(Hint)*T).*reshape(UI(:,end), [64,4]);
sqnorm(UIback)
max(abs(U(:,end) - UIback(:)))
UIback = exp(-1i*diag(Hint)*T).*reshape(UI(:,end), [64,4]);
max(abs(U(:,end) - UIback(:)))
est_errorsI
est_errors
busHamiltonian2qI
est_errorsI
matvecsI
matvecs
matvecsI
busHamiltonian2qI
est_errorsI
mniterI
matvecsI
matvecs
busHamiltonian2qI
[UIpar, mniterI, matvecsI, est_errorsI, historyI] = SemiGlobal1(GopI, Gdiff_opI, 0, [], [-50*1i, 50*1i], U2q0par(:), tgrid, 10e3, 9, 2, 1e-5, options, data);
busHamiltonian2q
tic, [Upar_ex, mniter_ex, matvecs_ex, est_errors_ex, history_ex] = SemiGlobal1(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0par(:), tgrid, 2.5e4, 9, 13, eps);toc
Upar_ex2qT = Upar_ex(:,end);
norm(U(:,end) - UIback(:))
norm(Upar_ex2qT - UIback(:))
norm(Upar_ex2qT - UIback(restore_vec))
norm(Upar_ex2qT - Upar(:,end))
est_errorsI
save Upar_ex_bus2q Upar_ex2qT
est_errors
norm(Upar_ex2qT - Upar(:,end))
norm(Upar_ex2qT - Upar(:,end))/norm(Upar_ex2qT)
norm(Upar_ex2qT)
sqnorm(Upar_ex2qT)
norm(Upar_ex2qT - UIback(restore_vec))/2
error = vecnorm(Upar - Upar_ex)/2;
figure
plot(tgrid, error)
figure
plot(tgrid(2:end), error(2:end) - error(1:end-1))
hold on
busHamiltonian2q
plot(tgrid(2:end), history.total_error)
size(history.total_error)
plot(tgrid(2:end), history.total_error(8:8:8000))
plot(tgrid(2:end), sum(sum(reshape(history.total_error, [8, 1000]))))
clf
plot(tgrid(2:end), error(2:end) - error(1:end-1))
hold on
plot(tgrid(2:end), sum(reshape(history.total_error, [8, 1000])))
hold on
plot(tgrid(2:end), cumsum(sum(reshape(history.total_error, [8, 1000]))))
est_errors
UIback_all = repmat(exp(-1i*diag(Hint)*T), [4,1]).*UI;
max(abs(UIback_all(:, end) - UIback(:)))
size(repmat(exp(-1i*diag(Hint)*T), [4,1]))
size(UIback_all)
size(UIback)
size(UIback(:))
UIback = exp(-1i*diag(Hint)*T).*reshape(UI(:,end), [64,4]);
max(abs(UIback_all(:, end) - UIback(:)))
norm(Upar_ex2qT - UIback_all(restore_vec, end))/2
est_errorsI
norm(Upar_ex2qT - Upar(:,end))/norm(Upar_ex2qT)
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop, Gdiff_op, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 6e3, 10, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_psi, Gdiff_op_psi, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 6e3, 10, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 6e3, 10, 1);
0.1^(1/18)
0.1^(1/14)
0.1^(1/12)
tic, [psip2, mniterp2, matvecsp2, est_errorsp2, historyp2] = SemiGlobal1(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), tgrid, 8e3, 9, 5, 1e-5);toc
f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts
(f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts)^(-1/14)
(10*f_scalar_error*max(abs(ev_domain))^Nt_ts/factorialNt_ts)^(-1/14)
factorial(14)/2*(4/530)^5
factorial(14)/2*(4/530)^5/530^9
(factorial(14)/2*(4/530)^5/530^9)^(1/14)
ans/T
T/(factorial(14)/2*(4/530)^5/530^9)^(1/14)
T/(0.1*factorial(14)/2*(4/530)^5/530^9)^(1/14)
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 4e3, 10, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 4.8e3, 10, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 5e3, 1, 1);
all_est_ers
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 6e3, 1, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 5.5e3, 1, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 5.8e3, 1, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 6e3, 1, 1);
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 6e3, 14, 1);
figure
plot(log10(allmv(1:10)), log10(aller(1:10)), '-o')
polyfit(log10(allmv(2:8)), log10(all_est_ers.conv_texp(2:8)), 1)
hold on
polyfit(log10(allmv(2:8)), log10(aller(2:8)), 1)
plot(log10(allmv(1:10)), log10(all_est_ers.conv(1:10)), '-o')
plot(log10(allmv(1:10)), log10(all_est_ers.texp_exact(1:10)), '-o')
plot(log10(allmv(1:10)), log10(all_est_ers.fm(1:10)), '-o')
[allNt97, allmv97, aller97, all_est_ers97] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 5e3, 10, 1);
[allNt97, allmv97, aller97, all_est_ers97] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 5e3, 1, 1);
[allNt97, allmv97, aller97, all_est_ers97] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 6e3, 1, 1);
[allNt97, allmv97, aller97, all_est_ers97] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 5.5e3, 1, 1);
[allNt97, allmv97, aller97, all_est_ers97] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 7, 5.5e3, 10, 1);
hold on
plot(log10(allmv(1:10)), log10(aller(1:10)), '-o')
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 6e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 6.5e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 7e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 8e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 6e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 5.5e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 5e3, 1, 1);
[allNt93, allmv93, aller93, all_est_ers93] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 3, 5e3, 14, 1);
plot(log10(allmv93(1:12)), log10(aller93(1:12)), '-o')
all_est_ers93
plot(log10(allmv93(1:12)), log10(all_est_ers.fm(1:12)), '-o')
plot(log10(allmv93(1:12)), log10(all_est_ers93.fm(1:12)), '-o')
[allNt94, allmv94, aller94, all_est_ers94] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 4, 5e3, 14, 1);
[allNt94, allmv94, aller94, all_est_ers94] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 4, 5e3, 12, 1);
plot(log10(allmv94), log10(aller94(1:12)), '-o')
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.5e4, 1);
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1e4, 1);
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.3e4, 1);
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.5e4, 15);
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.5e4, 25);
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.5e4, 20);
[allNtRK4, allmvRK4, allerRK4] =  error_decayRK4(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.5e4, 25);
figure
plot(log10(allmv(1:10)), log10(aller(1:10)), '-o')
hold on
plot(log10(allmv94), log10(aller94(1:12)), '-o')
plot(log10(allmv93(1:12)), log10(aller93(1:12)), '-o')
plot(log10(allmvRK4(1:22)), log10(allerRK4(1:22)), '-o')
polyfit(log10(allmvRK4(5:21)), log10(allerRK4(5:21)), 1)
allmvRK4(21)
allmvRK4(22)
allmvRK4(21)
allmv(5)
allmv(3)
allmvRK4(15)
[allNt95, allmv95, aller95, all_est_ers95] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], 9, 5, 6e3, 10, 1);
save bus2q_comparisons allNt95 allmv95 aller95 all_est_ers95 allNt94 allmv94 aller94 all_est_ers94 allNt93 allmv93 aller93 all_est_ers93 allNtRK4 allmvRK4 allerRK4
[allNtRK7, allmvRK7, allerRK7] =  error_decayRK7(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1e4, 1);
[allNtRK7, allmvRK7, allerRK7] =  error_decayRK7(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.3e4, 1);
[allNtRK7, allmvRK7, allerRK7] =  error_decayRK7(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.2e4, 1);
[allNtRK7, allmvRK7, allerRK7] =  error_decayRK7(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.1e4, 1);
[allNtRK7, allmvRK7, allerRK7] =  error_decayRK7(@(t, u) -1i*(H0odd*u - control(t)*(Ncoupler_even*u)), [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 1.2e4, 15);
polyfit(log10(allmvRK7(1:10)), log10(allerRK4(1:10)), 1)
polyfit(log10(allmvRK7(1:10)), log10(allerRK7(1:10)), 1)
plot(log10(allmvRK7(1:14)), log10(allerRK7(1:14)), '-o')
allmvRK7(5)
allmv95(2)
save bus2q_comparisons allNt95 allmv95 aller95 all_est_ers95 allNt94 allmv94 aller94 all_est_ers94 allNt93 allmv93 aller93 all_est_ers93 allNtRK4 allmvRK4 allerRK4 allNtRK7 allmvRK7 allerRK7
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 2e4, 1);
params(:, ti)
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 2e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 2.5e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 3e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 5e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 20, 5e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0even*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 20, 10e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 20, 3e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 3e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 2e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 5e3, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 7e3, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 7e3, 25);
plot(log10(allmvPWC), log10(allerPWC), '-o')
polyfit(log10(allmvPWC(1:4)), log10(allerPWC(1:4)), 1)
polyfit(log10(allmvPWC(6:end)), log10(allerPWC(6:end)), 1)
save bus2q_comparisons allNt95 allmv95 aller95 all_est_ers95 allNt94 allmv94 aller94 all_est_ers94 allNt93 allmv93 aller93 all_est_ers93 allNtRK4 allmvRK4 allerRK4 allNtRK7 allmvRK7 allerRK7 allNtPWC allmvPWC allerPWC
114.2*6000
114.2*6
70*1.4
figure
plot(tgrid, control(tgrid))
whos
GopI_odd = @(u, t, v) -1i*(H0odd*v - control(t)*Ncoupler_even*v) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*v) + g2*(bdag_a2dag_odd*v)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*v) + g2*(b_a2odd*v)));
GopI_odd = @(u, t, v) -1i*(H0odd*v - control(t)*Ncoupler_even*v + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*v) + g2*(bdag_a2dag_odd*v)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*v) + g2*(b_a2odd*v)));
Gdiff_opI = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_even*u1) + (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(g1*bdag_a1dag_odd*u1 + g2*bdag_a2dag_odd*u1) + (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(g1*b_a1odd*u1 + g2*b_a2odd*u1));
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 10e3, 1, 1);
Gdiff_opI_odd = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_even*u1) + (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(g1*bdag_a1dag_odd*u1 + g2*bdag_a2dag_odd*u1) + (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(g1*b_a1odd*u1 + g2*b_a2odd*u1));
Gdiff_opI = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_u*u1) +...
(exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(g1*bdag_a1dag_u*u1 + g2*bdag_a2dag_u*u1) +...
(exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(g1*b_a1u*u1 + g2*b_a2u*u1));
tic, [UIpar_ex, mniterI_ex, matvecsI_ex, est_errorsI_ex, historyI_ex] = SemiGlobal1(GopI, Gdiff_opI, 0, [], [-50*1i, 50*1i], U2q0par(:), tgrid, 3e4, 9, 13, eps);toc
tic, [UIpar_ex, mniterI_ex, matvecsI_ex, est_errorsI_ex, historyI_ex] = SemiGlobal1(GopI, Gdiff_opI, 0, [], [-50*1i, 50*1i], U2q0par(:), tgrid, 4e4, 9, 13, eps);toc
mniterI_ex
matvecsI_ex
UIpar_ex2qT = UIpar_ex(:,end);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 10e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 5, 10e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 11e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 12e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 17e3, 1, 1);
GopI_odd = @(u, t, v) -1i*(H0Iodd*v - control(t)*Ncoupler_even*v + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*v) + g2*(bdag_a2dag_odd*v)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*v) + g2*(b_a2odd*v)));
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 1.1e4, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 1e4, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 8e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 6e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 4e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 4.8e3, 1, 1);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 4.8e3, 1, 15);
[allNtI93, allmvI93, allerI93, all_est_ersI93] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 3, 4.8e3, 15, 1);
hold on
plot(log10(allmv95(1:10)), log10(aller95(1:10)), '-o')
figure
plot(log10(allmvI93), log10(allerI93), '-o')
polyfit(log10(allmv93(3:7)), log10(aller93(3:7)), 1)
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 1.5e4, 1);
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 1.2e4, 1);
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 1e4, 1);
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 5e3, 1);
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 3e3, 1);
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 4e3, 1);
[allNtIRK4, allmvIRK4, allerIRK4] =  error_decayRK4(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 3e3, 25);
hold on
plot(log10(allmvRK4(1:22)), log10(allerRK4(1:22)), '-o')
10^5.15
10^5.87
figure
plot(log10(allmvI93), log10(allerI93), '-o')
hold on
plot(log10(allmvIRK4(1:22)), log10(allerIRK4(1:22)), '-o')
plot(log10(allmvIRK4(1:23)), log10(allerIRK4(1:23)), '-o')
[allNtIPWC, allmvIPWC, allerIPWC] =  error_decayPWC(@(u, t) H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u)), [-50, 50], [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 10, 7e3, 1);
[allNtIPWC, allmvIPWC, allerIPWC] =  error_decayPWC(@(u, t) H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u)), [-50, 50], [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 10, 7e3, 25);
hold on
plot(log10(allmvPWC), log10(allerPWC), '-o')
plot(log10(allmvIPWC), log10(allerIPWC), '-o')
[allNtIRK7, allmvIRK7, allerIRK7] =  error_decayRK7(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 4e3, 1);
[allNtIRK7, allmvIRK7, allerIRK7] =  error_decayRK7(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 2e3, 1);
[allNtIRK7, allmvIRK7, allerIRK7] =  error_decayRK7(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 2e3, 20);
hold on
plot(log10(allmvRK7(1:14)), log10(allerRK7(1:14)), '-o')
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
[allNtI75, allmvI75, allerI75, all_est_ersI75] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 7, 5, 4.8e3, 1, 1);
[allNtI75, allmvI75, allerI75, all_est_ersI75] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 7, 5, 4.8e3, 15, 1);
tic, [psiI2_57, mniterI_57, matvecsI_57, est_errorsI_57, historyI_57] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), tgrid, 4.8e3, 7, 5, options_c, data_c);toc
options_c = options;
options_c.Niter = 1;
data_c = SGdata(options_c);
tic, [psiI2_57, mniterI_57, matvecsI_57, est_errorsI_57, historyI_57] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), tgrid, 4.8e3, 7, 5, [],  options_c, data_c);toc
norm(psiI2_57 - UIpar_ex2qT(33:64))
options_c.tol1st = eps;
tic, [psiI2_57, mniterI_57, matvecsI_57, est_errorsI_57, historyI_57] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 4.8e3, 7, 5, [],  options_c, data_c);toc
norm(psiI2_57(:, end) - UIpar_ex2qT(33:64))
est_errorsI_57
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 4.8e3, 7, 3, [],  options_c, data_c);toc
norm(psiI2_73(:, end) - UIpar_ex2qT(33:64))
est_errorsI_73
matvecsI_73
log10(43254)
tic, [psiI2_53, mniterI_53, matvecsI_53, est_errorsI_53, historyI_53] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 5, 3, [],  options_c, data_c);toc
norm(psiI2_53(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_53, mniterI_53, matvecsI_53, est_errorsI_53, historyI_53] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 5, 3, [],  options_c, data_c);toc
norm(psiI2_53(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_53, mniterI_53, matvecsI_53, est_errorsI_53, historyI_53] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 8e3, 5, 3, [],  options_c, data_c);toc
norm(psiI2_53(:, end) - UIpar_ex2qT(33:64))
est_errorsI_53
tic, [psiI2_53, mniterI_53, matvecsI_53, est_errorsI_53, historyI_53] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 8e3, 5, 2, [],  options_c, data_c);toc
tic, [psiI2_53, mniterI_53, matvecsI_53, est_errorsI_53, historyI_53] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 8e3, 5, 3, [],  options_c, data_c);toc
tic, [psiI2_52, mniterI_52, matvecsI_52, est_errorsI_52, historyI_52] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 8e3, 5, 2, [],  options_c, data_c);toc
est_errorsI_52
norm(psiI2_52(:, end) - UIpar_ex2qT(33:64))
10^(-3.986)
matvecsI_53
log10(56035)
allNtIRK7
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 7, 3, [],  options_c, data_c);toc
est_errorsI_73
norm(psiI2_73(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 8, 2, [],  options_c, data_c);toc
est_errorsI_73
norm(psiI2_73(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 5, 5, [],  options_c, data_c);toc
norm(psiI2_73(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 7, 3, [],  options_c, data_c);toc
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 9, 9, [],  options_c, data_c);toc
norm(psiI2_73(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_73, mniterI_73, matvecsI_73, est_errorsI_73, historyI_73] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 7, 3, [],  options_c, data_c);toc
tic, [psiI2_99, mniterI_99, matvecsI_99, est_errorsI_99, historyI_99] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 9, 9, [],  options_c, data_c);toc
est_errorsI_99
norm(psiI2_99(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_35, mniterI_35, matvecsI_35, est_errorsI_35, historyI_35] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 5e3, 3, 5, [],  options_c, data_c);toc
norm(psiI2_35(:, end) - UIpar_ex2qT(33:64))
est_errorsI_35
figure
plot(tgrid(2:end), historyI_73.reldif)
plot(tgrid(2:end), historyI_73.reldif(1:7:end))
plot(tgrid(2:end), historyI_73.reldif(1:6:end))
size(historyI_73.reldif)
plot(tgrid(2:end), historyI_73.reldif(5:5:end))
save bus2q_comparisons allNt95 allmv95 aller95 all_est_ers95 allNt94 allmv94 aller94 all_est_ers94 allNt93 allmv93 aller93 all_est_ers93 allNtRK4 allmvRK4 allerRK4 allNtRK7 allmvRK7 allerRK7 allNtPWC allmvPWC allerPWC allNtI93 allmvI93 allerI93 all_est_ersI93 allNtIRK4 allmvIRK4 allerIRK4 allNtIRK7 allmvIRK7 allerIRK7 allNtIPWC allmvIPWC allerIPWC
tic, [psiI2_99, mniterI_99, matvecsI_99, est_errorsI_99, historyI_99] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 9, 9, [],  options_c, data_c);toc
est_errorsI_99
norm(psiI2_99(:, end) - UIpar_ex2qT(33:64))
tic, [psi2_99, mniter_99, matvecs_99, est_errors_99, history_99] = SemiGlobal1(Gop_odd, Gdiff_op_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 9, 9, [],  options_c, data_c);toc
tic, [psi2_99, mniter_99, matvecs_99, est_errors_99, history_99] = SemiGlobal1(Gop_odd, Gdiff_op_even, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 9, 9, [],  options_c, data_c);toc
tic, [psi2_99, mniter_99, matvecs_99, est_errors_99, history_99] = SemiGlobal1(Gop_odd, Gdiff_op_even, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 9, 5, [],  options_c, data_c);toc
tic, [psi2_99, mniter_99, matvecs_99, est_errors_99, history_99] = SemiGlobal1(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), [0 T], 6e3, 9, 9, [],  options_c, data_c);toc
norm(psiI2_99(:, end) - Upar_ex2qT(33:64))
norm(psi2_99(:, end) - Upar_ex2qT(33:64))
figure
plot(tgrid(2:end), history_99.reldif(8:8:end))
size(history_99.reldif)
plot(tgrid(2:end), history_99.reldif(6:6:end))
figure
plot(tgrid(2:end), history_99.conv_error(6:6:end))
figure
plot(tgrid(2:end), history_99.conv_error(6:6:end)./history_99.reldif(6:6:end))
figure
plot(tgrid(2:end), historyI_99.reldif(6:6:end))
figure
plot(tgrid(2:end), historyI_99.conv_error(6:6:end))
1/0.8
10/13.5
figure
plot(tgrid(2:end), historyI_99.conv_error(6:6:end)./historyI_99.reldif(6:6:end))
mean(historyI_99.conv_error(6:6:end)./historyI_99.reldif(6:6:end))
mean(history_99.conv_error(6:6:end)./history_99.reldif(6:6:end))
p_max = 10; Mmax = 15;
all_gfuns = zeros(p_max + 1, Mmax - 1);
for Mi = 1:(Mmax - 1)
all_gfuns(:, Mi) = conv_gfuns(Mi + 1, p_max, 1);
end
all_gfuns
est_errorsI_99
tic, [psi2_89, mniter_89, matvecs_89, est_errors_89, history_89] = SemiGlobal1(Gop_odd, Gdiff_op_even, 0, [], [-528*1i, 1i], U2q0par(:, 2), [0 T], 6e3, 8, 9, [],  options_c, data_c);toc
est_errorsI_89
tic, [psiI2_89, mniterI_89, matvecsI_89, est_errorsI_89, historyI_89] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 9, 9, [],  options_c, data_c);toc
est_errorsI_89
norm(psiI2_89(:, end) - Upar_ex2qT(33:64))
norm(psiI2_89(:, end) - UIpar_ex2qT(33:64))
8/(4^8*(8^2-1)*(8-3))
sum(historyI_89.reldif)
MplusKmax = 30;
all_hfuns = zeros(p_max + 1, MplusKmax - 1);
for MplusKi = 1:(MplusKmax - 1)
all_hfuns(:, MplusKi) = conv_hfuns(MplusKi + 1, p_max, 1);
end
all_hfuns(:, 1:20)
all_guess_erI89 = historyI_89.texp_error_exact/A8*allg_funs(1,7);
A8 = 8/(4^8*(8^2-1)*(8-3));
all_guess_erI89 = historyI_89.texp_error_exact/A8*allg_funs(1,7);
all_guess_erI89 = historyI_89.texp_error_exact/A8*all_gfuns(1,7);
sum(all_guess_erI89)
sum(historyI_99.reldif)
sum(historyI_99.reldif) - sum(historyI_89.reldif)
tic, [psiI2_89, mniterI_89, matvecsI_89, est_errorsI_89, historyI_89] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 8, 9, [],  options_c, data_c);toc
est_errorsI_89
sum(historyI_89.reldif)
all_guess_erI89 = historyI_89.texp_error_exact/A8*allg_funs(1,7);
all_guess_erI89 = historyI_89.texp_error_exact/A8*all_gfuns(1,7);
sum(all_guess_erI89)
mean(historyI_99.conv_error(6:6:end)./historyI_99.reldif(6:6:end))
mean(historyI_89.conv_error(6:6:end)./historyI_89.reldif(6:6:end))
sum(all_guess_erI99)
tic, [psiI2_10_9, mniterI_10_9, matvecsI_10_9, est_errorsI_10_9, historyI_10_9] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 10, 9, [],  options_c, data_c);toc
est_errorsI_10_9
norm(psiI2_10_9(:, end) - UIpar_ex2qT(33:64))
mean(historyI_10_9.conv_error(6:6:end)./historyI_10_9.reldif(6:6:end))
sum(historyI_10_9.reldif)
sum(historyI_8_9.reldif)
sum(historyI_89.reldif)
sum(historyI_99.reldif)
all_guess_erI10_9 = historyI_10_9.texp_error_exact/A10*all_gfuns(1,7);
AM = M./(4.^M.*(M.^2-1).*(M-3));
M=1:15;
AM = M./(4.^M.*(M.^2-1).*(M-3));
A(8) - A8
A(8) - A7
A8
AM(8)-A8
AM
all_guess_erI10_9 = historyI_10_9.texp_error_exact/AM(10)*all_gfuns(1,9);
sum(all_guess_erI10_9)
tic, [psiI2_11_9, mniterI_11_9, matvecsI_11_9, est_errorsI_11_9, historyI_11_9] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 11, 9, [],  options_c, data_c);toc
est_errorsI_11_9
norm(psiI2_11_9(:, end) - UIpar_ex2qT(33:64))
mean(historyI_11_9.conv_error(6:6:end)./historyI_11_9.reldif(6:6:end))
sum(historyI_11_9.reldif)
sum(historyI_10_9.reldif)
tic, [psiI2_12_9, mniterI_12_9, matvecsI_12_9, est_errorsI_12_9, historyI_12_9] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 12, 9, [],  options_c, data_c);toc
est_errorsI_12_9
norm(psiI2_12_9(:, end) - UIpar_ex2qT(33:64))
T/6e3
T/6e3*50
pi/omega2
T/6e3/ans
tic, [psiI2_13_9, mniterI_13_9, matvecsI_13_9, est_errorsI_13_9, historyI_13_9] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 13, 9, [],  options_c, data_c);toc
est_errorsI_13_9
norm(psiI2_13_9(:, end) - UIpar_ex2qT(33:64))
tic, [psiI2_13_13, mniterI_13_13, matvecsI_13_13, est_errorsI_13_13, historyI_13_13] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 13, 13, [],  options_c, data_c);toc
est_errorsI_13_13
norm(psiI2_13_13(:, end) - UIpar_ex2qT(33:64))
clear psiI2_13_13 mniterI_13_13 matvecsI_13_13 est_errorsI_13_13 historyI_13_13
allNtIRK7
tic, [psiI2_12_2, mniterI_12_2, matvecsI_12_2, est_errorsI_12_2, historyI_12_2] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 6e3, 12, 2, [],  options_c, data_c);toc
est_errorsI_12_2
norm(psiI2_12_2(:, end) - UIpar_ex2qT(33:64))
allerIRK7
[allNtI12_2, allmvI12_2, allerI12_2, all_est_ersI12_2] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 12, 2, 4.8e3, 1, 1);
[allNtI12_2, allmvI12_2, allerI12_2, all_est_ersI12_2] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 12, 2, 5e3, 1, 1);
[allNtI12_2, allmvI12_2, allerI12_2, all_est_ersI12_2] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 12, 2, 6e3, 1, 1);
[allNtI12_2, allmvI12_2, allerI12_2, all_est_ersI12_2] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 12, 2, 5.5e3, 1, 1);
[allNtI12_2, allmvI12_2, allerI12_2, all_est_ersI12_2] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 12, 2, 6e3, 10, 1);
hold on
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
plot(log10(allmvI93), log10(allerI93), '-o')
[allNtI5_3_2i, allmvI5_3_2i, allerI5_3_2i, all_est_ersI5_3_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 5, 3, 4e3, 1, 2);
[allNtI5_3_2i, allmvI5_3_2i, allerI5_3_2i, all_est_ersI5_3_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 5, 3, 2e3, 1, 2);
[allNtI5_3_2i, allmvI5_3_2i, allerI5_3_2i, all_est_ersI5_3_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 5, 3, 2e3, 15, 2);
polyfit(log10(allmv5_3_2i(5:end)), log10(aller5_3_2i(5:end)), 1)
polyfit(log10(allmvI5_3_2i(5:end)), log10(allerI5_3_2i(5:end)), 1)
hold on
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
[allNtI7_3_2i, allmvI7_3_2i, allerI7_3_2i, all_est_ersI7_3_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 7, 3, 2e3, 15, 2);
hold on
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
polyfit(log10(allmvI7_3_2i(1:12)), log10(allerI7_3_2i(1:12)), 1)
hold on
plot(log10(allmvI93), log10(allerI93), '-o')
10^0.3
[allNtI9_2_2i, allmvI9_2_2i, allerI9_2_2i, all_est_ersI9_2_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 2, 2e3, 1, 2);
[allNtI9_2_2i, allmvI9_2_2i, allerI9_2_2i, all_est_ersI9_2_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 2, 3e3, 1, 2);
[allNtI9_2_2i, allmvI9_2_2i, allerI9_2_2i, all_est_ersI9_2_2i] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 9, 2, 3e3, 15, 2);
polyfit(log10(allmvI9_3_2i(3:7)), log10(allerI9_3_2i(3:7)), 1)
polyfit(log10(allmvI9_2_2i(3:7)), log10(allerI9_2_2i(3:7)), 1)
hold on
plot(log10(allmvI93), log10(allerI93), '-o')
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
10^0.2
10^0.13
allmvIRK7
[allmvIRK7 allerIRK7]
norm(UIpar(33:64) - UIpar_ex2qT(33:64))
norm(UIpar(:,33:64) - UIpar_ex2qT(33:64))
norm(UIpar(33:64, end) - UIpar_ex2qT(33:64))
matvecsI
busHamiltonian2qRWA
est_errorsRWA
busHamiltonian2qRWA
est_errorsRWA
infidelityRWA
infidelityI
sqnorm(UIpar(:,end) - URWApar(:,end))/4
reshape(URWA(:,end).*conj(URWA(:,end)), [64,4])
[reshape(UI(:,end).*conj(UI(:,end)), [64,4]), ans]
H0excitations = qubit_excitationsH(H0I, [4,4])
size(H0excitations)
Ncoupler_excitations = qubit_excitationsH(Ncoupler, [4,4])
eig(H0excitations)
max(control(tgrid))
min(control(tgrid))
eig(H0excitations-13*Ncoupler_excitations)
[URWA, mniterRWA, matvecsRWA, est_errorsRWA, historyRWA] = SemiGlobal1(@(u, t, v) -1i*(H0_excitations*v - control(t)*(Ncoupler_excitations*v)), @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_excitations*u1)), 0, [], [-30*1i, 20*1i], U2q0par(:), tgrid, 2e3, 9, 2, 1e-5, options, data);
[~, Nsingle, Ndouble] = qubit_excitationsH(H0I, [4,4])
u0RWA = zeros(13,1);
u0RWA([1, 2, 6, 8]) = 1
[excitations, indices] = excitation_comb(2, 3)
u0RWA = zeros(13,1);
u0RWA([1, 2, 6, 9]) = 1
[URWA, mniterRWA, matvecsRWA, est_errorsRWA, historyRWA] = SemiGlobal1(@(u, t, v) -1i*(H0_excitations*v - control(t)*(Ncoupler_excitations*v)), @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_excitations*u1)), 0, [], [-30*1i, 20*1i], u0RWA, tgrid, 2e3, 9, 2, 1e-5, options, data);
[URWA, mniterRWA, matvecsRWA, est_errorsRWA, historyRWA] = SemiGlobal1(@(u, t, v) -1i*(H0excitations*v - control(t)*(Ncoupler_excitations*v)), @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_excitations*u1)), 0, [], [-30*1i, 20*1i], u0RWA, tgrid, 2e3, 9, 2, 1e-5, options, data);
est_errorsRWA
singlei = excitation_kroni(1, [4,4])
doublei = excitation_kroni(2, [4,4])
URWAf = reshape(URWA(:, end), [64,4])
size(URWA(:, end))
URWAfull(restore_vec, :) = URWApar;
URWAfull_end = reshape(URWAfull(:, end), [64,4])
size(URWAfull(:, end))
URWAfull = zeros(256, 1001);
URWAfull(restore_vec, :) = URWApar;
URWAfull_end = reshape(URWAfull(:, end), [64,4])
URWAfull_end(singlei,2)
URWAfull_end(singlei,2)-URWA(2:4, end)
URWAfull_end(singlei,3)-URWA(5:7, end)
URWAfull_end(singlei,4)-URWA(8:13, end)
URWAfull_end(doublei,4)-URWA(8:13, end)
Usqrt_iswap = sparse(13);
Usqrt_iswap(indices2q, indices2q) = Usqrt_iswap_2q;
indices2q = [1, 2, Nsingle + 3, 2*Nsingle + 3]
URWAsqrt_iswap = sparse(Ntotal);
URWAsqrt_iswap(indices2q, indices2q) = Usqrt_iswap_2q
Ntotal = 1 + 2*Nsingle + Ndouble
URWAsqrt_iswap = sparse(Ntotal);
URWAsqrt_iswap(indices2q, indices2q) = Usqrt_iswap_2q
URWAold = URWA;
fidelityRWAold = fidelityRWA;
busHamiltonian2qRWA
size(URWA(:,end))
size(URWAsqrt_iswap(:))
URWAsqrt_iswap = sparse(Ntotal, 1);
URWAsqrt_iswap([1, 2:3, Nsingle + 2:3, 2*Nsingle + 3]) = Usqrt_iswap_2q([1, 6:7, 10:11, 16])
URWAsqrt_iswap([1, 2:3, Nsingle + (2:3), 2*Nsingle + 3]) = Usqrt_iswap_2q([1, 6:7, 10:11, 16])
bus
busHamiltonian2qRWA
max(max(abs(URWA-URWAold)))
fidelityRWA-fidelityRWAold
fidelityRWA
H0excitations(2:4, 2:4)
full(H0excitations(2:4, 2:4))
eig(ans)
H0single = H0excitations(2:4, 2:4);
Ncoupler_single = Ncoupler_excitations(2:4, 2:4);
Gop_single = @(u, t, v) -1i*(H0single*v - control(t)*(Ncoupler_single*v));
Gdiff_op_single = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_single*u1));
tic, [psiRWA2, mniterRWA2, matvecsRWA2, est_errorsRWA2, historyRWA2] = SemiGlobal1(GopRWA, Gdiff_opRWA, 0, [], [-16*1i, 9*1i], u0RWA(2:4), tgrid, 2e3, 9, 2, [], options_c, data_c); toc
tic, [psiRWA2, mniterRWA2, matvecsRWA2, est_errorsRWA2, historyRWA2] = SemiGlobal1(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), tgrid, 2e3, 9, 2, [], options_c, data_c); toc
est_errorsRWA2
matvecsRWA2
matvecsRWA
mniterRWA
mniterRWA2
[allNtRWA92, allmvRWA92, allerRWA92, all_est_ersRWA92] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 9, 2, 1e3, 1, 1);
tic,[URWAex, mniterRWAex, matvecsRWAex, est_errorsRWAex, historyRWAex] = SemiGlobal1(GopRWA, Gdiff_opRWA, 0, [], [-30*1i, 20*1i], u0RWA, tgrid, 15e3, 9, 13, eps, options, data);toc
est_errorsRWAex
size(UIpar_ex2qT)
URWA_ex2qT = URWAex(:,end);
figure
plot(tgrid, URWAex.*conj(URWAex))
[allNtRWA92, allmvRWA92, allerRWA92, all_est_ersRWA92] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 9, 2, 1e3, 1, 1);
[allNtRWA92, allmvRWA92, allerRWA92, all_est_ersRWA92] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 9, 2, 1.5e3, 1, 1);
[allNtRWA92, allmvRWA92, allerRWA92, all_est_ersRWA92] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 9, 2, 1.2e3, 1, 1);
[allNtRWA92, allmvRWA92, allerRWA92, all_est_ersRWA92] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 9, 2, 1.1e3, 1, 1);
[allNtRWA92, allmvRWA92, allerRWA92, all_est_ersRWA92] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 9, 2, 1.1e3, 15, 1);
[allNtRWA_RK7, allmvRWA_RK7, allerRWA_RK7] =  error_decayRK7(@(t, u) -1i*(H0single*u - control(t)*(Ncoupler_single*u)), [0 T], u0RWA(2:4), URWA_ex2qT(2:4), 1e3, 1);
[allNtRWA_RK7, allmvRWA_RK7, allerRWA_RK7] =  error_decayRK7(@(t, u) -1i*(H0single*u - control(t)*(Ncoupler_single*u)), [0 T], u0RWA(2:4), URWA_ex2qT(2:4), 0.5e3, 1);
[allNtRWA_RK7, allmvRWA_RK7, allerRWA_RK7] =  error_decayRK7(@(t, u) -1i*(H0single*u - control(t)*(Ncoupler_single*u)), [0 T], u0RWA(2:4), URWA_ex2qT(2:4), 0.5e3, 15);
[allNtRWA_RK7, allmvRWA_RK7, allerRWA_RK7] =  error_decayRK7(@(t, u) -1i*(H0single*u - control(t)*(Ncoupler_single*u)), [0 T], u0RWA(2:4), URWA_ex2qT(2:4), 0.5e3, 20);
hold on
plot(log10(allmvRWA_RK7(1:16)), log10(allerRWA_RK7(1:16)), '-o')
all_est_ersRWA92
all_est_ersRWA92.fm
all_est_ersRWA92.conv
[all_est_ersRWA92.fm.' all_est_ersRWA92.conv]
[allNtRWA52, allmvRWA52, allerRWA52, all_est_ersRWA52] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 2, 1.1e3, 1, 1);
[allNtRWA52, allmvRWA52, allerRWA52, all_est_ersRWA52] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 2, 1e3, 1, 1);
[allNtRWA52, allmvRWA52, allerRWA52, all_est_ersRWA52] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 2, 0.5e3, 1, 1);
[allNtRWA52, allmvRWA52, allerRWA52, all_est_ersRWA52] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 2, 0.7e3, 1, 1);
[allNtRWA52, allmvRWA52, allerRWA52, all_est_ersRWA52] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 2, 0.7e3, 15, 1);
polyfit(log10(allmvRWA52(2:end)), log10(allerRWA52(2:end)), 1)
polyfit(log10(allmvRWA52(3:end)), log10(allerRWA52(3:end)), 1)
hold on
plot(log10(allmvRWA52), log10(all_est_ersRWA52.texp_exact), '-o')
hold on
plot(log10(allmvRWA52), log10(all_est_ersRWA52.conv), '-o')
plot(log10(allmvRWA52), log10(all_est_ersRWA52.fm), '-o')
figure
plot(log10(allmvRWA52), log10(allerRWA52), '-o')
hold on
plot(log10(allmvRWA_RK7(1:16)), log10(allerRWA_RK7(1:16)), '-o')
[allNtRWA53, allmvRWA53, allerRWA53, all_est_ersRWA53] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 3, 0.7e3, 17, 1);
polyfit(log10(allmvRWA53(3:end)), log10(allerRWA53(3:end)), 1)
hold on
plot(log10(allmvRWA_RK7(1:16)), log10(allerRWA_RK7(1:16)), '-o')
H0single
full(H0excitations(2:4, 2:4))
eig(ans)
omega1-omega2
[allNtRWA55, allmvRWA55, allerRWA55, all_est_ersRWA55] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 5, 5, 0.7e3, 17, 1);
polyfit(log10(allmvRWA55(3:end)), log10(allerRWA55(3:end)), 1)
hold on
plot(log10(allmvRWA_RK7(1:16)), log10(allerRWA_RK7(1:16)), '-o')
[allNtRWA73, allmvRWA73, allerRWA73, all_est_ersRWA73] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 7, 3, 0.7e3, 1, 1);
[allNtRWA73, allmvRWA73, allerRWA73, all_est_ersRWA73] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 7, 3, 1e3, 1, 1);
[allNtRWA73, allmvRWA73, allerRWA73, all_est_ersRWA73] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 7, 3, 1e3, 17, 1);
polyfit(log10(allmvRWA73(2:12)), log10(allerRWA73(2:12)), 1)
hold on
plot(log10(allmvRWA_RK7(1:16)), log10(allerRWA_RK7(1:16)), '-o')
plot(log10(allmvRWA92), log10(allerRWA92), '-o')
figure
plot(log10(allmvRWA73), log10(allerRWA73), '-o')
plot(log10(allmvRWA73(1:13)), log10(allerRWA73(1:13)), '-o')
hold on
plot(log10(allmvRWA73(1:13)), log10(all_est_ersRWA73.fm(1:13)), '-o')
plot(log10(allmvRWA73(1:13)), log10(all_est_ersRWA73.conv(1:13)), '-o')
T/(2*pi)
T*omega_phi
10^0.8
T/6e3*2*omega_phi/(2*pi)
T/6e3*2*omega1/(2*pi)
PWCersI2_99 = PWCerrors(historyI2_99.U(:, 1:8:end), GopIpar, ev_domain, T, 20);
whos
size(historyI_99.U)
PWCersI2_99 = PWCerrors(historyI_99.U(:, 1:8:end), GopI_odd, [-50*1i, 50*1i], T, 20);
GopI_odd
3+5+4+1+5+ (17+10+49+23+32+21)/60
(3+5+4+1+5+ (17+10+49+23+32+21)/60)*20
(750*3+1050*4+850)/8
(750*3+1050*4+850)/8 + (3+5+4+1+5+ (17+10+49+23+32+21)/60)*20
[allNtI82, allmvI82, allerI82, all_est_ersI82] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 8, 2, 3e3, 15, 1);
[allNtI82, allmvI82, allerI82, all_est_ersI82] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 8, 2, 4e3, 1, 1);
[allNtI82, allmvI82, allerI82, all_est_ersI82] = error_decaySG2(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), UIpar_ex2qT(33:64), [0, T], 8, 2, 4e3, 15, 1);
polyfit(log10(allmvI82(3:11)), log10(allerI82(3:11)), 1)
hold on
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
[allNtIRK7b, allmvIRK7b, allerIRK7b] =  error_decayRK7(@(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))), [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 4e3, 13);
hold on
plot(log10(allmvIRK7b), log10(allerIRK7b), '-o')
[allNtRWA82, allmvRWA82, allerRWA82, all_est_ersRWA82] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 8, 2, 1e3, 1, 1);
[allNtRWA82, allmvRWA82, allerRWA82, all_est_ersRWA82] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 8, 2, 1e3, 17, 1);
polyfit(log10(allmvRWA82(2:10)), log10(allerRWA82(2:10)), 1)
[allNtRWA_RK7b, allmvRWA_RK7b, allerRWA_RK7b] =  error_decayRK7(@(t, u) -1i*(H0single*u - control(t)*(Ncoupler_single*u)), [0 T], u0RWA(2:4), URWA_ex2qT(2:4), 1e3, 17);
hold on
plot(log10(allmvRWA_RK7b), log10(allerRWA_RK7b), '-o')
allNtIRK7b
all_ers = PWCerrors(historyI_ex.U(33:64, 1:10:end), HopRK, T/8000:T/4000:(T-T/8000), [-9, 16], [0 T], 20);
HopRK = @(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u)))
all_ers4e3 = PWCerrors(historyI_ex.U(33:64, 1:10:end), HopRK, T/8000:T/4000:(T-T/8000), [-9, 16], [0 T], 20);
GopI_RK = @(t, u) -1i*(H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u)))
clear HopRK
all_ers4e3 = PWCerrors(historyI_ex.U(33:64, 1:10:end), HopRK, T/8000:T/4000:(T-T/8000), [-50, 50], [0 T], 20);
HopIPWC = @(u, t) H0Iodd*u - control(t)*(Ncoupler_even*u) + exp(1i*2*omega1*t)*(g1*(bdag_a1dag_odd*u) + g2*(bdag_a2dag_odd*u)) + exp(-1i*2*omega1*t)*(g1*(b_a1odd*u) + g2*(b_a2odd*u))
all_ers4e3 = PWCerrors(historyI_ex.U(33:64, 1:10:end), HopIPWC, T/8000:T/4000:(T-T/8000), [-50, 50], [0 T], 20);
all_ers4e3 = PWCerrors(historyI_ex.U(33:64, 1:80:end), HopIPWC, T/8000:T/4000:(T-T/8000), [-50, 50], [0 T], 20);
figure
plot((1:4e3)*T/4000,  all_ers4e3)
sum(all_ers4e3)
allerI82
all_gfuns
all_hfuns
[allNtRWA83, allmvRWA83, allerRWA83, all_est_ersRWA83] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 8, 3, 1e3, 13, 1);
polyfit(log10(allmvRWA83(2:10)), log10(allerRWA83(2:10)), 1)
[allNtRWA84, allmvRWA84, allerRWA84, all_est_ersRWA84] = error_decaySG2(Gop_single, Gdiff_op_single, 0, [], [-16*1i, 9*1i], u0RWA(2:4), URWA_ex2qT(2:4), [0, T], 8, 4, 1e3, 13, 1);
polyfit(log10(allmvRWA84(2:10)), log10(allerRWA84(2:10)), 1)
hold on
plot(log10(allmvRWA82), log10(allerRWA_RK7b), '-o')
plot(log10(allmvRWA82), log10(allerRWA82), '-o')
plot(log10(allmvRWA82(1:13)), log10(allerRWA82(1:13)), '-o')
plot(log10(allmvRWA84), log10(allerRWA84), '-o')
figure
plot(log10(allmvRWA82(1:13)), log10(allerRWA82(1:13)), '-o')
hold on
plot(log10(allmvRWA82(1:13)), log10(all_est_ersRWA82.conv(1:13)), '-o')
polyfit(log10(allmvRWA82(2:10)), log10(all_est_ersRWA82.conv(2:10)), 1)
polyfit(log10(allmvRWA82(3:10)), log10(all_est_ersRWA82.conv(3:10)), 1)
polyfit(log10(allmvRWA82(3:8)), log10(all_est_ersRWA82.conv(3:8)), 1)
all_est_ersRWA82.fm
all_est_ersRWA83.fm
all_est_ersRWA83.texp_exact
all_est_ersRWA82.texp_exact
T/1e3/(2*pi)
T/1e3/(2*pi)*omega_phi
whos
size(historyI.U)
size(historyIex.U)
size(historyI_ex.U)
norm(Dchebb(historyI_ex.U(:, 160001:160009)))
D2M = D2chebbMat(9, T/4e4);
vecnorm((D2M*Dchebb(historyI_ex.U(33:64, 160001:160009)).').')
vecnorm((D2M*(historyI_ex.U(33:64, 160001:160009)).').')
size(D2M)
D2M = D2chebbMat(8, T/4e4);
vecnorm((D2M*(historyI_ex.U(33:64, 160001:160009)).').')
DM = DchebbMat(8, T/4e4);
vecnorm((DM*(historyI_ex.U(33:64, 160001:160009)).').')
T/1e3
T/1e3*vecnorm((D2M*(historyI_ex.U(33:64, 160001:160009)).').')
vecnorm((DM^4*DM*(historyI_ex.U(33:64, 160001:160009)).').')
vecnorm((DM^4*DM*(historyI_ex.U(33:64, 160001:160009)).').')*T/4e3
vecnorm((DM^4*(historyI_ex.U(33:64, 160001:160009)).').')
vecnorm((DM^4*(historyI_ex.U(33:64, 240001:240009)).').')
vecnorm((DM^4*DM*(historyI_ex.U(33:64, 240001:240009)).').')*T/4e3
all_ersS4e3 = PWCerrors(history_ex.U(33:64, 1:80:end), HopPWC, T/8000:T/4000:(T-T/8000), [-1, 530], [0 T], 20);
HopPWC = @(u, t) H0odd*u - control(t)*(Ncoupler_even*u)
all_ersS4e3 = PWCerrors(history_ex.U(33:64, 1:80:end), HopPWC, T/8000:T/4000:(T-T/8000), [-1, 530], [0 T], 20);
figure
plot((1:4e3)*T/4000,  all_ersS4e3)
all_ersS4e3 = PWCerrors(history_ex.U(33:64, 1:40:end), HopPWC, T/5000:T/5000:(T-T/5000), [-1, 530], [0 T], 20);
all_ersS4e3 = PWCerrors(history_ex.U(33:64, 1:40:end), HopPWC, T/10000:T/5000:(T-T/10000), [-1, 530], [0 T], 20);
plot((1:4e3)*T/4000,  all_ersS4e3)
plot((1:5e3)*T/5000,  all_ersS4e3)
tic, [psiI2_82_4e3, mniterI82_4e3, matvecsI82_4e3, est_errorsI82_4e3, historyI82_4e3] = SemiGlobal1(GopI_odd, Gdiff_opI_odd, 0, [], [-50*1i, 50*1i], U2q0par(:, 2), [0 T], 4e3, 8, 2, [],  options_c, data_c);toc
est_errorsI82_4e3
figure
plot((1:4e3)*T/4000,  all_ersS4e3)
plot((1:4e3)*T/4000,  all_ers4e3)
hold on
plot((1:4e3)*T/4000,  historyI82_4e3.reldif)
all_ers0th4e3 = vecnorm(historyI82_4e3.U(8:7:end)-historyI82_4e3.U(1:7:(end-7)))./vecnorm(historyI82_4e3.U(8:7:end));
size(all_ers0th4e3)
all_ers0th4e3 = vecnorm(historyI82_4e3.U(:, 8:7:end)-historyI82_4e3.U(:, 1:7:(end-7)))./vecnorm(historyI82_4e3.U(:,8:7:end));
size(all_ers0th4e3)
plot((1:4e3)*T/4000,  all_ers0th4e3)
2*omega1
figure
plot(log10(allmv93), log10(aller93), '-o')
plot(log10(allmv93(1:12)), log10(aller93(1:12)), '-o')
hold on
plot(log10(allmv95(1:10)), log10(aller95(1:10)), '-o')
plot(log10(allmvRK4(1:22)), log10(allerRK4(1:22)), '-o')
plot(log10(allmvRK7(1:14)), log10(allerRK7(1:14)), '-o')
plot(log10(allmvPWC), log10(allerPWC), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
figure
plot(log10(allmvI82), log10(allerI82))
figure
plot(log10(allmvI82), log10(allerI82), '-o')
hold on
plot(log10(allmvI92), log10(allerI92), '-o')
plot(log10(allmvI93), log10(allerI93), '-o')
plot(log10(allmvI), log10(allerI), '-o')
whos
plot(log10(allmvIRK4(1:23)), log10(allerIRK4(1:23)), '-o')
plot(log10(allmvIRK7b), log10(allerIRK7b), '-o')
doc save
save bus2q_comparisons allNtRKI7b allmvRKI7b allerRKI7b allNtI82 allmvI82 allerI82 all_est_ersI82 -append
save bus2q_comparisons allNtIRK7b allmvIRK7b allerIRK7b allNtI82 allmvI82 allerI82 all_est_ersI82 -append
plot(log10(allmvIPWC), log10(allerIPWC), '-o')
plot(log10(allmvIRK7(1:16)), log10(allerIRK7(1:16)), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
UIRK2 = RK2(GopRK, [0, T], U2q0par(:), 1e4);
UIRK2 = RK2(GopIRK, [0, T], U2q0par(:), 1e4);
UIRK2 = RK2(GopI_RK, [0, T], U2q0par(:), 1e4);
norm(UIRK2 - UIpar_ex2qT)
UIRK2 = RK2(GopI_RK, [0, T], U2q0par(:), 3e4);
norm(UIRK2 - UIpar_ex2qT)
clear UIRK2
psi2IRK2 = RK2(GopI_RK, [0, T], U2q0par(33:64), 1e4);
clear psi2IRK2
psiI2RK2 = RK2(GopI_RK, [0, T], U2q0par(33:64), 1e4);
norm(psiI2RK2(:, end) - UIpar_ex2qT(33:64))
psiI2RK2 = RK2(GopI_RK, [0, T], U2q0par(33:64), 3e4);
norm(psiI2RK2(:, end) - UIpar_ex2qT(33:64))
psiI2RK2 = RK2(GopI_RK, [0, T], U2q0par(:, 2), 1e4);
norm(psiI2RK2(:, end) - UIpar_ex2qT(33:64))
norm(psiI2RK2(:, end))
psiI2RK2(:, end)
psiI2RK2 = RK2(GopI_RK, [0, T], U2q0par(:, 2), T/1e4);
norm(psiI2RK2(:, end))
psiI2RK2 = RK2(GopI_RK, [0, T], U2q0par(:, 2), T/3e4);
norm(psiI2RK2(:, end))
norm(psiI2RK2(:, end) - UIpar_ex2qT(33:64))
psifI2RK2 = RK2uf(GopI_RK, [0, T], U2q0par(:, 2), T/3e4);
max(max(abs(psiI2RK2(:, end) - psifI2RK2)))
[allNtIRK2, allmvIRK2, allerIRK2] =  error_decayRK7(GopI_RK, [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), , 13);
[allNtIRK2, allmvIRK2, allerIRK2] =  error_decayRK7(GopI_RK, [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 7e3, 25);
[allNtIRK2, allmvIRK2, allerIRK2] =  error_decayRK2(GopI_RK, [0 T], U2q0par(:, 2), UIpar_ex2qT(33:64), 7e3, 25);
hold on
plot(log10(allmvIPWC), log10(allerIPWC), '-o')
figure
plot(log10(allNtIRK2), log10(allerIRK2), '-o')
hold on
plot(log10(allNtIPWC), log10(allerIPWC), '-o')
allNtIPWC./allNtIRK2
allerIPWC./allerIRK2
H0double = H0excitations(7:13, 7:13);
eig(H0double)
Gop_double = @(u, t, v) -1i*(H0double*v - control(t)*(Ncoupler_double*v));
Ncoupler_double = Ncoupler_excitations(7:13, 7:13);
Gop_double = @(u, t, v) -1i*(H0double*v - control(t)*(Ncoupler_double*v));
Gdiff_op_double = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_double*u1));
tic, [psiRWA4_82, mniterRWA4_82, matvecsRWA4_82, est_errorsRWA4_82, historyRWA4_82] = SemiGlobal1(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(7:13), tgrid, 2e3, 8, 2, [], options_c, data_c); toc
mniterRWA4_82
est_errorsRWA4_82
est_errorsRWA2
norm(psiRWA4_82(:,end) -URWA(7:13, end))
norm(psiRWA2(:,end) -URWA(2:4, end))
tic, [psiRWA4_92, mniterRWA4_92, matvecsRWA4_92, est_errorsRWA4_92, historyRWA4_92] = SemiGlobal1(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(7:13), tgrid, 2e3, 8, 2, 1e-5, options, data); toc
est_errorsRWA4_82
est_errorsRWA4_92
mniterRWA4_92
norm(psiRWA4_92(:,end) -URWA(7:13, end))
Ncoupler_double = Ncoupler_excitations(7:13, 7:13);
Gop_double = @(u, t, v) -1i*(H0double*v - control(t)*(Ncoupler_double*v));
Gdiff_op_double = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_double*u1));
u0RWA
H0double = H0excitations(8:13, 8:13);
Ncoupler_double = Ncoupler_excitations(8:13, 8:13);
Gop_double = @(u, t, v) -1i*(H0double*v - control(t)*(Ncoupler_double*v));
Gdiff_op_double = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_double*u1));
tic, [psiRWA4_92, mniterRWA4_92, matvecsRWA4_92, est_errorsRWA4_92, historyRWA4_92] = SemiGlobal1(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(7:13), tgrid, 2e3, 9, 2, 1e-5, options, data); toc
tic, [psiRWA4_92, mniterRWA4_92, matvecsRWA4_92, est_errorsRWA4_92, historyRWA4_92] = SemiGlobal1(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(8:13), tgrid, 2e3, 9, 2, 1e-5, options, data); toc
est_errorsRWA4_92
mniterRWA4_92
norm(psiRWA4_92(:,end) -URWA(8:13, end))
tic, [psiRWA4_82, mniterRWA4_82, matvecsRWA4_82, est_errorsRWA4_82, historyRWA4_82] = SemiGlobal1(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(8:13), tgrid, 2e3, 8, 2, [], options_c, data_c); toc
est_errorsRWA4_82
norm(psiRWA4_82(:,end) -URWA(8:13, end))
[allNtRWA4_83, allmvRWA4_83, allerRWA4_83, all_est_ersRWA4_83] = error_decaySG2(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(8:13), URWA_ex2qT(8:13), [0, T], 8, 3, 1e3, 1, 1);
[allNtRWA4_83, allmvRWA4_83, allerRWA4_83, all_est_ersRWA4_83] = error_decaySG2(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(8:13), URWA_ex2qT(8:13), [0, T], 8, 3, 1.5e3, 1, 1);
all_est_ersRWA4_83
[allNtRWA4_82, allmvRWA4_82, allerRWA4_82, all_est_ersRWA4_82] = error_decaySG2(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(8:13), URWA_ex2qT(8:13), [0, T], 8, 2, 1.5e3, 1, 1);
all_est_ersRWA4_82
[allNtRWA4_82, allmvRWA4_82, allerRWA4_82, all_est_ersRWA4_82] = error_decaySG2(Gop_double, Gdiff_op_double, 0, [], [-30*1i, 20*1i], u0RWA(8:13), URWA_ex2qT(8:13), [0, T], 8, 2, 1.5e3, 13, 1);
hold on
plot(log10(allmvRWA82(1:13)), log10(all_est_ersRWA82.conv(1:13)), '-o')
[allNtRWA4_RK7, allmvRWA4_RK7, allerRWA4_RK7] =  error_decayRK7(@(t, u) -1i*(H0double*u - control(t)*(Ncoupler_double*u)), [0 T], u0RWA(8:13), URWA_ex2qT(8:13), 1e3, 17);
hold on
plot(log10(allmvRWA_RK7b), log10(allerRWA_RK7b), '-o')
plot(log10(allmvRWA4_RK7b), log10(allerRWA4_RK7b), '-o')
plot(log10(allmvRWA4_RK7), log10(allerRWA4_RK7), '-o')
plot(log10(allmvRWA_RK7b), log10(allerRWA_RK7b), '-o')
historyRWA4_82
historyRWA
tic, [psiRWA4_82ar, mniterRWA4_82ar, matvecsRWA4_82ar, est_errorsRWA4_82ar, historyRWA4_82ar] = SemiGlobal1(Gop_double, Gdiff_op_double, 0, [], [], u0RWA(8:13), tgrid, 2e3, 8, 2, [], options_c, data_c); toc
est_errorsRWA4_82ar
est_errorsRWA4_82
historyRWA4_82ar
figure
plot((1:2000)*T/2e3, historyRWA4_82ar.stab_factor)
min(historyRWA4_82ar.stab_factor)
max(historyRWA4_82ar.stab_factor)
whos
save Upar_ex_bus2q URWA_ex2qT -append
save Upar_ex_bus2q UIpar_ex2qT -append
save bus2q_comparisons allNtRWA_RK7b allmvRWA_RK7b allerRWA_RK7b allNtRWA82 allmvRWA82 allerRWA82 all_est_ersRWA82 -append
% The energy units are g1=2*pi*100MHz
omegac0 = 7445/100;
omega1 = 5889.9/100;
omega2 = 5031.1/100;
omega3 = 4350/100;
alphac = 230/100;
alpha1 = 324/100;
alpha2 = 235/100;
alpha3 = 1;
g1 = 1;
g2 = 71.4/100;
g3 = 85/100;
n = (0:3).';
N_HO = spdiags(n, 0, 4, 4);
Manhar = spdiags([0; 0; -1; -3], 0, 4, 4);
Mcoup = spdiags(sqrt([(n + 1), n]), [-1, 1], 4, 4);
I4 = speye(4);
% The order of the sub-systems in the direct product:
% Tunable bus coupler, transmon 3, transmon 1, transmon 2
H03q = multi_kron({omegac0*N_HO + alphac*Manhar, I4, I4, I4}) + multi_kron({I4, I4, omega1*N_HO + alpha1*Manhar, I4})...
+ multi_kron({I4, I4, I4, omega2*N_HO + alpha2*Manhar}) + multi_kron({I4, omega3*N_HO + alpha3*Manhar, I4, I4})...
+ g1*multi_kron({Mcoup, I4, Mcoup, I4}) + g2*multi_kron({Mcoup, I4, I4, Mcoup}) + g3*multi_kron({Mcoup, Mcoup, I4, I4});
size(H03q)
eig(H03q)
eig(H0)
Ncoupler3q = multi_kron({N_HO, I4, I4, I4});
max(control(tgrid))
min(control(tgrid))
eig(H03q-13*Ncoupler3q)
busHamiltonian3q
size(Upar3q)
mniter3q
est_errors3q
matvecs3q
0.1^(1/13)
8e3/ans
busHamiltonian3q
matvecs3q
mniter3q
est_errors3q
infidelity3q
infidelity2q
infidelity
busHamiltonian2q
busHamiltonian3q
matvecs3q
mniter3q
est_errors3q
busHamiltonian3q
est_errors3q
norm(Upar3q(:, end))
history3q.stab_factor
mniter3q
busHamiltonian3q
mniter3q
est_errors3q
busHamiltonian3q
matvecs3q
est_errors3q
tic, [Upar3q_ex, mniter3q_ex, matvecs3q_ex, est_errors3q_ex, history_ex] = SemiGlobal1(Gop3q, Gdiff_op3q, 0, [], [-655*1i, 1i], U2q0par3q(:), tgrid, 3e4, 9, 13, eps, options, data); toc
Upar_ex3qT = Upar3q_ex(:, end);
reshape(U3q(:,end).*conj(U3q(:,end)), [256,4])
norm(Upar_ex3qT - Upar3q(:,end))/norm(Upar_ex3qT)
save bus2q_comparisons Gop_odd Gop_even GopI_odd GopI_even -append
save bus2q_comparisons Gop_odd Gop_even GopI_odd -append
save bus2q_comparisons GopI_RK -append
save bus2q_comparisons Gdiff_op_even Gdiff_opI_even -append
save bus2q_comparisons Gdiff_op_even Gdiff_opI_odd -append
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff3q_op_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 6e3, 1, 1);
Gop3q_odd = @(u, t, v) -1i*(H03q_odd*v - control(t)*Ncoupler_even*v);
Gdiff_op_odd = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler_even*u1);
Gdiff_op3q_odd = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler3q_even*u1);
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 6e3, 1, 1);
Gop3q_odd = @(u, t, v) -1i*(H03q_odd*v - control(t)*Ncoupler3q_even*v);
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 6e3, 1, 1);
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 8e3, 1, 1);
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 7e3, 1, 1);
all_est_ers3q95
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 7.5e3, 1, 1);
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 7.2e3, 1, 1);
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 7.5e3, 1, 1);
all_est_ers3q95
[allNt3q95, allmv3q95, aller3q95, all_est_ers3q95] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 9, 5, 7.5e3, 12, 1);
polyfit(log10(allmv3q95(2:7)), log10(aller3q95(2:7)), 1)
[allNt3qRK7, allmv3qRK7, aller3qRK7] =  error_decayRK7(@(t, u) -1i*(H03qodd*u - control(t)*(Ncoupler3q_even*u)), [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 1.2e4, 1);
[allNt3qRK7, allmv3qRK7, aller3qRK7] =  error_decayRK7(@(t, u) -1i*(H03q_odd*u - control(t)*(Ncoupler3q_even*u)), [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 1.2e4, 1);
[allNt3qRK7, allmv3qRK7, aller3qRK7] =  error_decayRK7(@(t, u) -1i*(H03q_odd*u - control(t)*(Ncoupler3q_even*u)), [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 1.3e4, 1);
[allNt3qRK7, allmv3qRK7, aller3qRK7] =  error_decayRK7(@(t, u) -1i*(H03q_odd*u - control(t)*(Ncoupler3q_even*u)), [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 1.4e4, 1);
[allNt3qRK7, allmv3qRK7, aller3qRK7] =  error_decayRK7(@(t, u) -1i*(H03q_odd*u - control(t)*(Ncoupler3q_even*u)), [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 1.4e4, 12);
hold on
plot(log10(allmv3qRK7), log10(aller3qRK7), '-o')
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 5, 7.5e3, 1, 1);
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 5, 7e3, 1, 1);
all_est_ers3q85
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 5, 6e3, 1, 1);
all_est_ers3q85
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 6, 6e3, 1, 1);
all_est_ers3q85
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 6, 5e3, 1, 1);
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 6, 5.5e3, 1, 1);
all_est_ers3q85
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(Gop3q_odd, Gdiff_op3q_odd, 0, [], [-655*1i, 1i], U2q0par3q(:, 2), Upar_ex3qT(129:256), [0, T], 8, 6, 5.5e3, 13, 1);
allNt3q86 = allNt3q85; allmv3q86=allmv3q85; aller3q86=aller3q85; all_est_ers3q86=all_est_ers3q85;
clear allNt3q85 allmv3q85 aller3q85 all_est_ers3q85
polyfit(log10(allmv3q86(2:8)), log10(aller3q86(2:8)), 1)
plot(log10(allmv3q86), log10(aller3q86), '-o')
Gop3qRK = @(t, u) -1i*(H03q_odd*u - control(t)*(Ncoupler3q_even*u))
hold on
plot(log10(allmv3q86), log10(all_est_ers3q86.conv), '-o')
plot(log10(allmv3q86), log10(all_est_ers3q86.fm), '-o')
plot(log10(allmv3q86), log10(all_est_ers3q86.texp_error), '-o')
plot(log10(allmv3q86), log10(all_est_ers3q86.texp_exact), '-o')
figure
plot(log10(allmv3q95), log10(aller3q95), '-o')
hold on
plot(log10(allmv3q95), log10(all_est_ers3q95.conv), '-o')
plot(log10(allmv3q95), log10(all_est_ers3q95.fm), '-o')
plot(log10(allmv3q95), log10(all_est_ers3q95.texp_exact), '-o')
[allNt3qRK4, allmv3qRK4, aller3qRK4] =  error_decayRK4(Gop3qRK, [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 2e4, 1);
[allNt3qRK4, allmv3qRK4, aller3qRK4] =  error_decayRK4(Gop3qRK, [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 1.4e4, 1);
[allNt3qRK4, allmv3qRK4, aller3qRK4] =  error_decayRK4(Gop3qRK, [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 2e4, 22);
plot(log10(allmv3qRK4), log10(aller3qRK4), '-o')
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 10, 7e3, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 10, 8e3, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 10, 1e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 15, 1e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 20, 1e4, 1);
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 20, 1e4, 15);
[allNtPWC10, allmvPWC10, allerPWC10] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 10, 1e4, 15);
10^0.3
[allNtPWC10, allmvPWC10, allerPWC10] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u), [-655*1i, 1i], [0 T], U2q0par3q(:, 2), Upar_ex3qT(129:256), 10, 2e4, 15);
hold on
plot(log10(allmvPWC), log10(allerPWC), '-o')
allNt3qPWC20 = allNtPWC; allmv3qPWC20 = allmvPWC; aller3qPWC20 = allerPWC;
allNt3qPWC10 = allNtPWC10; allmv3qPWC10 = allmvPWC10; aller3qPWC10 = allerPWC10;
clear allNtPWC10 allmvPWC10 allerPWC10
[allNtPWC, allmvPWC, allerPWC] =  error_decayPWC(@(u, t) H0odd*u - control(t)*(Ncoupler_even*u), [-1, 528], [0 T], U2q0par(:, 2), Upar_ex2qT(33:64), 10, 7e3, 25);
plot(log10(allmv3qPWC10), log10(aller3qPWC10), '-o')
plot(log10(allmv3qPWC20), log10(aller3qPWC20), '-o')
save bus3q_comparisons allNt3q95 allmv3q95 aller3q95 all_est_ers3q95 allNt3q86 allmv3q86 aller3q86 all_est_ers3q86 allNt3qRK7 allmv3qRK7 aller3qRK7 allNt3qRK4 allmv3qRK4 aller3qRK4 allNt3qPWC10 allmv3qPWC10 aller3qPWC10 allNt3qPWC20 allmv3qPWC20 aller3qPWC20
Dc0 = omegac0 - omega1;
D2 = omega2 - omega1;
D3 = omega3 - omega1;
H0I3q = multi_kron({Dc0*N_HO + alphac*Manhar, I4, I4, I4}) + multi_kron({I4, I4, alpha1*Manhar, I4})...
+ multi_kron({I4, I4, I4, D2*N_HO + alpha2*Manhar}) + multi_kron({I4, D3*N_HO + alpha3*Manhar, I4, I4})...
+ g1*multi_kron({Mcoup, I4, Mcoup, I4}) + g2*multi_kron({Mcoup, I4, I4, Mcoup}) + g3*multi_kron({Mcoup, Mcoup, I4, I4});
eig(H0I3q)
eig(H0I3q-13*Ncoupler3q)
busHamiltonian3qI
infidelityI3q
fidelityI3q
fidelity3q
norm(UI3q(:, end))
2-norm(UI3q(:, end))
est_errorsI3q
mniterI3q
infidelityI
infidelity
busHamiltonian3qI
infidelityI3q
infidelity3q
norm(UI3q(:, end))
reshape(U3q(:,end).*conj(U3q(:,end)), [256,4])
reshape(UI3q(:,end).*conj(UI3q(:,end)), [256,4])
busHamiltonian3qI
infidelity3q
infidelityI3q
max(abs(U(:,end) - UIback(:)))
Hint3q = multi_kron({I4, I4, omega1*N_HO, I4}) + multi_kron({omega1*N_HO, I4, I4, I4}) + multi_kron({I4, I4, I4, omega1*N_HO}) + multi_kron({I4, omega1*N_HO, I4, I4});
UIback3q = exp(1i*diag(Hint)*T).*reshape(UI(:,end), [64,4]);
UIback3q = exp(-1i*diag(Hint)*T).*reshape(UI(:,end), [64,4]);
UIback3q = exp(-1i*diag(Hint)*T).*reshape(UI(:,end), [256,4]);
UIback3q = exp(-1i*diag(Hint)*T).*reshape(UI3q(:,end), [256,4]);
UIback3q = exp(-1i*diag(Hint3q)*T).*reshape(UI3q(:,end), [256,4]);
UIback = exp(-1i*diag(Hint)*T).*reshape(UI(:,end), [64,4]);
UIback3q = exp(-1i*diag(Hint3q)*T).*reshape(UI3q(:,end), [256,4]);
max(abs(U3q(:,end) - UIback3q(:)))
est_errorsI3q
mniterI3q
busHamiltonian3qI
mniterI3q
est_errorsI3q
matvecs3q
matvecsI3q
UIback3q = exp(-1i*diag(Hint3q)*T).*reshape(UI3q(:,end), [256,4]);
max(abs(U3q(:,end) - UIback3q(:)))
tic, [Upar3q_ex, mniter3q_ex, matvecs3q_ex, est_errors3q_ex, history_ex] = SemiGlobal1(GopI3q, Gdiff_opI3q, 0, [], [-42*1i, 95*1i], U2q0par3q(:), tgrid, 4e4, 9, 13, eps, options, data); toc
est_errors3q_ex
UIpar3q_ex = Upar3q_ex; mniterI3q_ex = mniter3q_ex, matvecsI3q_ex=matvecs3q_ex; est_errorsI3q_ex=est_errors3q_ex; historyI3q_ex = history_ex;
tic, [Upar3q_ex, mniter3q_ex, matvecs3q_ex, est_errors3q_ex, history_ex] = SemiGlobal1(Gop3q, Gdiff_op3q, 0, [], [-655*1i, 1i], U2q0par3q(:), tgrid, 3e4, 9, 13, eps, options, data); toc
tic, [Upar3q_ex, mniter3q_ex, matvecs3q_ex, est_errors3q_ex, history3q_ex] = SemiGlobal1(Gop3q, Gdiff_op3q, 0, [], [-655*1i, 1i], U2q0par3q(:), tgrid, 3e4, 9, 13, eps, options, data); toc
UIpar3q_old = UIpar3q;
atanh(1)
atanh(-1)
atanh(77)
busHamiltonian3qI
max(max(abs(UIpar3q_old - UIpar3q)))
UIpar_ex3qT = UIpar3q_ex(:, end);
norm(UIpar3q(:, end) - UIpar_ex3qT)
est_errorsI3q
GopI3q_odd = @(u, t, v) -1i*(H03q_odd*v - control(t)*Ncoupler3q_even*v + exp(1i*2*omega1*t)*(Mdouble_creation_odd*v) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*v));
Mdouble_creation_odd = g1*bdag_a1dag_3q_odd + g2*bdag_a2dag_3q_odd + g3*bdag_a3dag_3q_odd
Mdouble_annihilation_odd = g1*b_a1_3q_odd + g2*b_a2_3q_odd + g3*b_a3_3q_odd;
GopI3q_odd = @(u, t, v) -1i*(H03q_odd*v - control(t)*Ncoupler3q_even*v + exp(1i*2*omega1*t)*(Mdouble_creation_odd*v) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*v));
Gdiff_op
Gdiff_opI3q_odd = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler3q_odd*u1) + (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(Mdouble_creation_odd*u1) + (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(Mdouble_annihilation_odd*u1));
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 2, 1e-5, options, data);toc
Gdiff_opI3q_odd = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler3q_even*u1) + (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(Mdouble_creation_odd*u1) + (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(Mdouble_annihilation_odd*u1));
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 2, 1e-5, options, data);toc
mniter2I3q
norm(psi2I3q(:,end)
norm(psi2I3q(:,end))
GopI3q_odd = @(u, t, v) -1i*(H03q_odd*v - control(t)*Ncoupler3q_even*v + exp(1i*2*omega1*t)*(Mdouble_creation_odd*v) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*v));
Gdiff_opI3q_odd = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler3q_even*u1) + (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(Mdouble_creation_odd*u1) + (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(Mdouble_annihilation_odd*u1));
Mdouble_creation_odd - Mdouble_creation_odd'
Mdouble_creation - Mdouble_creation'
Mdouble_creation - Mdouble_annihilation'
Mdouble_creation_odd - Mdouble_annihilation_odd'
Ncoupler3q_even
U2q0par3q(:, 2)
U2q0par3q(:)
U2q0par3q
U2q0par3q(:,2)
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 2e4, 9, 5, 1e-5, options, data);toc
norm(psi2I3q(:,end))
est_errors2I3q
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 5, 1e-5, options, data);toc
est_errors2I3q
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 2, 1e-5, options, data);toc
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 3, 1e-5, options, data);toc
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 4, 1e-5, options, data);toc
est_errors2I3q
norm(psi2I3q(:,end))
mniter2I3q
norm(psi2I3q(:,end) - UIpar3q(129:256,end))
Mdouble_creation_odd - Mdouble_creation(129:256, 129:256)
Mdouble_annihilation_odd - Mdouble_annihilation(129:256, 129:256)
GopI3q_odd = @(u, t, v) -1i*(H03q_odd*v - control(t)*(Ncoupler3q_even*v) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*v) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*v));
Gdiff_opI3q_odd = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler3q_even*u1) + (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(Mdouble_creation_odd*u1) + (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(Mdouble_annihilation_odd*u1));
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 4, 1e-5, options, data);toc
norm(psi2I3q(:,end) - UIpar3q(129:256,end))
norm(psi2I3q(:,end))
1-norm(psi2I3q(:,end))
U2q0par3q
GopI3q_odd = @(u, t, v) -1i*(H0I3q_odd*v - control(t)*(Ncoupler3q_even*v) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*v) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*v));
tic, [psi2I3q, mniter2I3q, matvecs2I3q, est_errors2I3q, history2I3q] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.2e4, 9, 2, 1e-5, options, data);toc
est_errors2I3q
norm(psi2I3q(:,end) - UIpar3q(129:256,end))
[allNt3q85, allmv3q85, aller3q85, all_est_ers3q85] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 1e4, 1, 1);
[allNtI3q82, allmvI3q82, allerI3q82, all_est_ersI3q82] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 1e4, 1, 1);
all_est_ersI3q82
[allNtI3q82, allmvI3q82, allerI3q82, all_est_ersI3q82] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 0.9e4, 1, 1);
all_est_ersI3q82
[allNtI3q82, allmvI3q82, allerI3q82, all_est_ersI3q82] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 8e3, 1, 1);
all_est_ersI3q82
[allNtI3q82, allmvI3q82, allerI3q82, all_est_ersI3q82] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 7e3, 1, 1);
all_est_ersI3q82
[allNtI3q82, allmvI3q82, allerI3q82, all_est_ersI3q82] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 5e3, 1, 1);
all_est_ersI3q82
[allNtI3q82, allmvI3q82, allerI3q82, all_est_ersI3q82] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 8, 2, 5e3, 15, 1);
tic
[Upar3q, mniter3q, matvecs3q, est_errors3q, history3q] = SemiGlobal1(Gop3q, Gdiff_op3q, 0, [], [-655*1i, 1i], U2q0par3q(:), tgrid, 1e4, 9, 5, 1e-5, options, data);
toc
polyfit(log10(allmvI3q82(3:9)), log10(aller3qI3q82(3:9)), 1)
polyfit(log10(allmvI3q82(3:9)), log10(allerI3q82(3:9)), 1)
plot(log10(allmvI3q82), log10(all_est_ersI3q82.conv), '-o')
plot(log10(allmvI3q82), log10(allersI3q82), '-o')
plot(log10(allmvI3q82), log10(allerI3q82), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
hold on
plot(log10(allmvI3q82), log10(all_est_ersI3q82.conv), '-o')
polyfit(log10(allmvI3q82(3:9)), log10(all_est_ersI3q82.conv(3:9)), 1)
hold on
plot(log10(allmvI3q82), log10(all_est_ersI3q82.fm), '-o')
plot(log10(allmvI3q82), log10(all_est_ersI3q82.texp_exact), '-o')
figure
plot(log10(allmvI82), log10(allerI82), '-o')
hold on
plot(log10(allmvI82), log10(all_est_ersI82.conv), '-o')
hold on
polyfit(log10(allmvI82(4:10)), log10(all_est_ersI82.conv(4:10)), 1)
polyfit(log10(allmvI82(4:10)), log10(allerI82.conv(4:10)), 1)
polyfit(log10(allmvI82(4:10)), log10(allerI82(4:10)), 1)
polyfit(log10(allmvI82(4:10)), log10(allerI82.conv_cheb(4:10)), 1)
polyfit(log10(allmvI82(4:10)), log10(all_est_ersI82.conv_cheb(4:10)), 1)
polyfit(log10(allmvI82(4:10)), log10(all_est_ersI82.conv_texp(4:10)), 1)
polyfit(log10(allmvI82(4:10)), log10(all_est_ersI82.conv_fm(4:10)), 1)
plot(log10(allmvI82), log10(all_est_ersI82.conv_texp), '-o')
GopI3qRK = @(t, u) -1i*(H0I3q_odd*u - control(t)*(Ncoupler3q_even*u) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*u) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*u))
[allNtI3qRK7, allmvI3qRK7, allerI3qRK7] =  error_decayRK7(GopI3qRK, [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 2.5e3, 16);
plot(log10(allmvI3qRK7), log10(allerI3qRK7), '-o')
polyfit(log10(allmvI3qRK7(1:13)), log10(aller3qRK7(1:13)), 1)
polyfit(log10(allmvI3qRK7(1:13)), log10(allerI3qRK7(1:13)), 1)
plot(log10(allmvI3q82*1.5), log10(allerI3q82), '-o')
allNtI3qRK7
allNtI3qRK7(8)
allNtI3q82(5)
tic, [psi2I3q82, mniter2I3q82, matvecs2I3q82, est_errors2I3q82, history2I3q82] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c, data_c);toc
mniter2I3q82
tic, [psi2I3q82, mniter2I3q82, matvecs2I3q82, est_errors2I3q82, history2I3q82] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c, data_c);toc
tic, psi2I3qRK7 = RK7uf(GopI3qRK, [0 T], U2q0par3q(:, 2), T/1.25e4);toc
[psi2I3q82, mniter2I3q82, matvecs2I3q82, est_errors2I3q82, history2I3q82] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c, data_c);
norm(UIpar_ex3qT(129:256)-psi2I3qRK7)
psi2I3q82 = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c, data_c);
psi2I3qRK7 = RK7uf(GopI3qRK, [0 T], U2q0par3q(:, 2), T/1.25e4);
psi2I3q82 = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c, data_c);
options_c
options_c2 = options_c
options_c2.conv_er_cheb = false
tic, [psi2I3q82, mniter2I3q82, matvecs2I3q82, est_errors2I3q82, history2I3q82] = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c2, data_c);toc
tic, psi2I3q82 = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c2, data_c);toc
psi2I3q82 = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c2, data_c);
psi2I3q82 = SemiGlobal1(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), [0 T], 1.25e4, 8, 2, [], options_c, data_c);
[allNtI3q92, allmvI3q92, allerI3q92, all_est_ersI3q92] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 9, 2, 5e3, 1, 1);
[allNtI3q92, allmvI3q92, allerI3q92, all_est_ersI3q92] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 9, 2, 6e3, 1, 1);
[allNtI3q92, allmvI3q92, allerI3q92, all_est_ersI3q92] = error_decaySG2(GopI3q_odd, Gdiff_opI3q_odd, 0, [], [-42*1i, 95*1i], U2q0par3q(:, 2), UIpar_ex3qT(129:256), [0, T], 9, 2, 6e3, 13, 1);
plot(log10(allmvI3q92(1:9)), log10(all_est_ersI3q92(1:9)), '-o')
plot(log10(allmvI3q92(1:9)), log10(allerI3q92(1:9)), '-o')
polyfit(log10(allmvI3q92(2:7)), log10(allerI3q92(2:7)), 1)
[allNtI3qRK4, allmvI3qRK4, allerI3qRK4] =  error_decayRK4(GopI3qRK, [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 2.5e3, 1);
[allNtI3qRK4, allmvI3qRK4, allerI3qRK4] =  error_decayRK4(GopI3qRK, [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 2.5e3, 20);
plot(log10(allmvI3qRK4), log10(allerI3qRK4), '-o')
[allNtI3qPWC10, allmvI3qPWC10, allerI3qPWC10] =  error_decayPWC(@(u, t) H03q_odd*u - control(t)*(Ncoupler3q_even*u) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*u) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*u), [-42*1i, 95*1i], [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 10, 2e4, 1);
[allNtI3qPWC10, allmvI3qPWC10, allerI3qPWC10] =  error_decayPWC(@(u, t) H0I3q_odd*u - control(t)*(Ncoupler3q_even*u) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*u) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*u), [-42*1i, 95*1i], [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 10, 2e4, 1);
[allNtI3qPWC10, allmvI3qPWC10, allerI3qPWC10] =  error_decayPWC(@(u, t) H0I3q_odd*u - control(t)*(Ncoupler3q_even*u) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*u) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*u), [-42*1i, 95*1i], [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 10, 1e4, 1);
[allNtI3qPWC10, allmvI3qPWC10, allerI3qPWC10] =  error_decayPWC(@(u, t) H0I3q_odd*u - control(t)*(Ncoupler3q_even*u) + exp(1i*2*omega1*t)*(Mdouble_creation_odd*u) + exp(-1i*2*omega1*t)*(Mdouble_annihilation_odd*u), [-42*1i, 95*1i], [0 T], U2q0par3q(:, 2), UIpar_ex3qT(129:256), 10, 1e4, 20);
hold on
plot(log10(allmvI3qPWC10), log10(allerI3qPWC10), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')