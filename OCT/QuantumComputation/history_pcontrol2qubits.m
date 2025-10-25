alpha = abs(Enl(3))
normal_factor = 4/alpha
Enlg = Enl*normal_factor
thresholdg = 5*normal_factor
E7g = E7*normal_factor
Delta = 13;
[E7g; []]
[E7g.'; []]
isa([], 'function_handle')
save pcontrol2qubits
M = rand(1e3, 1e4);
for j = 1:1e3, L=M(:,end:-1:1); end
tic, for j = 1:100, L=M(:,end:-1:1); end, toc
tic, for j = 1:100, Q=fliplr(M); end, toc
max(max(abs(Q-L))))
max(max(abs(Q-L)))
clear M L Q
clear all
test_soft_int
whos
allJmax_overlap4cInt
test_soft_int
temp = allJmax_overlap4cInt;
test_soft_int
allJmax_overlap4cInt(1:2)
temp(1:2)-ans
temp(1:2
temp(1:2)
temp
clear all
integw_t = chebgridw(4, 7, 1);
integw_t
sum(integw_t(2:7))
clear all
load transmonH25
load pcontrol2qubits
whos
pE7
doc speye
speye(3)
I7s = speye(7)
aq
adagq
whos
Enlg
Manhar = spdiags(Enlg, 0, 7, 7)
H0 = multi_kron{I7s, Manhar, I7s} + multi_kron{I7s, I7s, Manhar}
H0 = multi_kron({I7s, Manhar, I7s}) + multi_kron({I7s, I7s, Manhar})
7^3
I5s = speye(5)
H05 = multi_kron({I5s, Manhar(1:5, 1:5), I5s}) + multi_kron({I5s, I5s, Manhar(1:5,1:5)})
p5 = 1i*(aqdag(1:5,1:5) - aq(1:5,1:5))
p5 = 1i*(adagq(1:5,1:5) - aq(1:5,1:5))
Hd1 = multi_kron({I5s, p5, I5s})
Hd2 = multi_kron({I5s, I5s, p5})
a
a7s = spdiags(sqrt((0:6).'), 1, 7, 7)
adag3s = spdiags(sqrt((1:7).'), -1, 7, 7)
clear adag3s
adag7s = spdiags(sqrt((1:7).'), -1, 7, 7)
Hcplus = multi_kron({a7s(1:5, 1:5), adagq(1:5, 1:5), I5s}) + multi_kron({a7s(1:5, 1:5), I5s, adagq(1:5, 1:5)})
Hcminus = Hcplus';
Hcminus - multi_kron({adag(1:,1:5), aq(1:5,1:5), I5s}) - multi_kron({adag(1:,1:5), I5s, aq(1:5,1:5)})
Hcminus - multi_kron({adag(1:5,1:5), aq(1:5,1:5), I5s}) - multi_kron({adag(1:5,1:5), I5s, aq(1:5,1:5)})
Hcminus - multi_kron({adag7(1:,1:5), aq(1:5,1:5), I5s}) - multi_kron({adag7(1:,1:5), I5s, aq(1:5,1:5)})
Hcminus - multi_kron({adag7(1:5,1:5), aq(1:5,1:5), I5s}) - multi_kron({adag7(1:5,1:5), I5s, aq(1:5,1:5)})
Hcminus - multi_kron({adag7s(1:5,1:5), aq(1:5,1:5), I5s}) - multi_kron({adag7s(1:5,1:5), I5s, aq(1:5,1:5)})
H05u = kron(speye(4), H05)
I4s = speye(4)
Hd1s = kron(I4s,Hd1)
Hd2s = kron(I4s,Hd2);
clear Hd1s Hd2s
Hd1u = kron(I4s,Hd1);
Hd2u = kron(I4s,Hd2);
Hcplusu = kron(I4s,Hcplus);
Hcminusu = kron(I4s,Hcminus);
generateH1modep
gs5 = zeros(5,1);
%-- 11/02/2021 18:25 --%
doc sparse
Utarget = [multi_kron({gs5, gs5, gs5}); multi_kron({gs5, gs5, fundamental1}); multi_kron({gs5, fundamental1, gs5});...
multi_kron({gs5, fundamental1, fundamental1})]
generateH1modep
Utarget
generateH1modep
Utarget
generateH1modep
U0
penalforbv
multi_kron({ones5, ones5, phi4});
multi_kron({ones5, ones5, phi4})
ones5
ones5(:) = 1
generateH1modep
ones5
penalforbv
penalforbv.*ones(125,1)
generateH1modep
penalforbv
(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}))
(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}))+0
ans.*U0
double(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}))
whos
nnz(penalforbv)
45*8
45*16
generateH1modep
penalforbv
save pcontrol2qubits
13*pi
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0{:}, Utarget{:}, (1:4)*125, fJpsit, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, Edomain, fguess, @(w) 0.5*(1-tanh(2*(w-13))), options, 13*pi, 1e3, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:4)*125, fJpsit, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, Edomain, fguess, @(w) 0.5*(1-tanh(2*(w-13))), options, 13*pi, 1e3, 7, 7, 1e-6);
figure
x=0:0.1:20;
plot(x, 0.5*(1-tanh(2*(x-13)))
plot(x, 0.5*(1-tanh(2*(x-13))))
options = optionsOCqn(1e-4, 1e4);
filterE = @(w) 0.5*(1-tanh(2*(w-13)))
fJforb = get_fJforb_penal(penalforbv)
eig(H05 + Hcplus + Hcminus)
eig(H05)
min(ans)
eig(H05 - Hcplus - Hcminus)
eig(H05 + Hcplus + Hcminus)
options.f_max_alpha = get_f_max_alphaOCf(5, 13*pi/1e3, 13*pi, filterE);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJpsit, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.5*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.5*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w).*0.5*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 7, 7, 1e-6);
size(penalforbv)
penalforbv = kron(ones(4,1),penalforbv);
penalforbv = double(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}));
penalforbvu = kron(ones(4,1),penalforbv);
fJforb = get_fJforb_penal(penalforbvu);
size(penalforbvu))
size(penalforbvu)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 7, 7, 1e-6);
size((ifft_v_ext(1:N, :).*exp(1i*wdt*normalig.')))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 9, 7, 1e-6);
options.f_max_alpha = get_f_max_alphaOCf(5, 13*pi/1e3, 13*pi, filterE);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 9, 7, 1e-6);
options = optionsOCqn(1e-4, 1e4);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 9, 7, 1e-6);
options.f_max_alpha = get_f_max_alphaOCf(5, 13*pi/1e3, 13*pi, filterE);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 9, 9, 1e-6);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/1e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1e3, 9, 9, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1200, 7, 7, 1e-6
options = optionsOCqn(1e-4, 1e4);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1200, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1500, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 2000, 7, 7, 1e-6);
sqnorm(allpsi(:,end))
4-ans
sqnorm(allchi(:,1))
sqnorm(allchi(:,end))
sqnorm(allchi(:,3000))
chiihterm(:,1)
max(max(abs(chiihterm)))
penalforbvu
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 2000, 7, 7, 1e-6);
max(max(abs(all_penal_psi)))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1000, 7, 7, 1e-6);
fJforb = get_fJforb_penal(penalforbvu);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1000, 7, 7, 1e-6);
fJforb = get_fJforb_penal(penalforbvu);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1000, 7, 7, 1e-6);
mallniterc2mb1
mallniterc
clear invHess
whos
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1200, 7, 7, 1e-6);
mallniterc
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), filterE, options, 13*pi, 1200, 7, 7, 1e-6);
fieldwnzM
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e-3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1200, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.1*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1200, 7, 7, 1e-6);
Jterms
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1200, 7, 7, 1e-6);
isreal(Jenergy_fun)
isreal(Jenergy)
isreal(Jenergy_fun)
isreal(Jenergy)
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 2000, 7, 7, 1e-6);
isreal(Jenergy)
isreal(Jenergy_fun)
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 2000, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1000, 7, 7, 1e-6);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/1e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1000, 7, 7, 1e-6);
fJforb = get_fJforb_penal(1e-3*penalforbvu);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1000, 7, 7, 1e-6);
Jterms
mallniterc
fJforb = get_fJforb_penal(1e-3*penalforbvu);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/1.2e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1200, 7, 7, 1e-6);
max(abs(allpsi(5,:)))
max(abs(Hparams(1,:))
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e2*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1200, 7, 7, 1e-6
fJforb = get_fJforb_penal(penalforbvu);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e2*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 1200, 7, 7, 1e-6);
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
max(abs(allpsi(5,:)))
max(abs(Hparams(2,:)))
max(abs(Hparams(1,:)))
%-- 14/02/2021 7:38 --%
load pcontrol2qubits
whos
x=0:0.1:20;
figure
plot(x, 0.5*(1-tanh(2*(x-13))))
plot(x, 0.5*(1-tanh(2*(x-13))).*heaviside(15-x))
edit heaviside
penalforbv = double(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}));
penalforbvu = kron(ones(4,1),penalforbv);
fJforb = get_fJforb_penal(penalforbvu);
options = optionsOCqn(1e-4, 1e4);
filterE = @(w) 0.5*(1-tanh(2*(w-13))).*heaviside(15-x);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.05*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))), @(w) 1e3*0.5.*(1-tanh(2*(w-13))), options, 13*pi, 2000, 7, 7, 1e-6)
fJforb = get_fJforb_penal(10*penalforbvu);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.05*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.05*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
filterE
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.05*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
filterE = @(w) 0.5*(1-tanh(2*(w-13))).*heaviside(15-w);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.05*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
mallniterc
max(abs(fieldt), 2)
max(abs(fieldt), [], 2)
max(abs(psi(5,:)))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w/13)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
max(abs(allpsi(5,:)))
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e2*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
Jterms
[fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, [], Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e2*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
max(abs(allpsi(5,:)))
max(abs(Hparams(2,:)))
max(abs(Hparams(1,:)))
max(abs(allpsi(5,:)))
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
[fieldt1, fieldw1, psi1, relE1, conv1, niter1, mallniterc1, Jterms1, maxgrad1, alpha1, invHess1] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, [], Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw, @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
max(abs(Hparams(2,:)))
max(abs(Hparams(1,:)))
max(abs(allpsi(5,:)))
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
[~, Jforb] = fJforb(allpsi, Nt_ts, integw, penalforbv)
penalforbv = double(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}));
penalforbvu = kron(ones(4,1),penalforbv);
figure
6.5/0.075
plot(0:13*pi/2e3:13*pi, fieldt)
figure
plot(0:1/13:2e3/13, fieldw)
plot(0:1/13:200/13, fieldw(1:201))
plot(0:1/13:200/13, fieldw(:,1:201))
projh5 = multi_kron({phi4, ones5, ones5});
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, [], Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw1, @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6
options.invHess0 = invHess1;
options = optionsOCqn(1e-4, 1e4);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
options1=options;
options1.invHess0 = invHess1;
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, [], Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw1, @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options1, 13*pi, 2000, 7, 7, 1e-6
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, [], Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw1, @(w) 1e3*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options1, 13*pi, 2000, 7, 7, 1e-6);
max(abs(Hparams(2,:)))
max(abs(Hparams(1,:)))
max(abs(allpsi(5,:)))
figure
whos
save pcontrol2qubits
figure
plot(0:1/13:200/13, fieldw2(:,1:201))
figure
plot(0:13*pi/2e3:13*pi, fieldt)
integw_t = chebgridw(2e3, 7, 13*pi/2e3);
[~, Jforb2] = fJforb(allpsi2, 7, integw_t, penalforbvu)
[~, Jforb2] = fJforb(psi2, 7, ones(1, 2001)*13*pi/2e3, penalforbvu)
fJforb_penal =fJforb
clear fJforb
[~, Jforb2] = fJforb(psi2, 7, ones(1, 2001)*13*pi/2e3, penalforbvu)
exval_projforb = real(sum(conj(psi2).*penalforbvu.*psi2));
figure
plot(0:13*pi/2e3:13*pi, exval_projforb)
plot(0:13*pi/2e3:13*pi, exval_projforb/4)
exval_projh4 = real(sum(conj(psi2).*projh4u.*psi2));
projh4 = zeros(125,1);
projh4(101:125) = 1;
projh4u=kron(ones(4,1), projh4);
save pcontrol2qubits
exval_projh4 = real(sum(conj(psi2).*projh4u.*psi2));
figure
plot(0:13*pi/2e3:13*pi, exval_projh4/4)
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb_penal, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw2, @(w) 1e4*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
fJforb1 = get_fJforb_penal(penalforbvu);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb1, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw2, @(w) 1e4*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
fJforb01 = get_fJforb_penal(0.1*penalforbvu);
[fieldt2, fieldw2, psi2, relE2, conv2, niter2, mallniterc2, Jterms2, maxgrad2, alpha2, invHess2] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb01, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw2, @(w) 1e4*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
max(abs(allpsi(5,:)))
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
max(abs(Hparams(1,:)))
max(abs(allpsi(5,:)))
1-Jmax
max(abs(allpsi(5,:)))
max(abs(Hparams(1,:)))
max(abs(Hparams(2,:)))
save pcontrol2qubits
Jterms2
1-Jterms.max
1-Jterms.Jmax
1-Jterms2.Jmax
exval_projh4 = real(sum(conj(psi2).*projh4u.*psi2));
figure
plot(0:13*pi/2e3:13*pi, exval_projh4/4)
exval_projforb = real(sum(conj(psi2).*penalforbvu.*psi2));
plot(0:13*pi/2e3:13*pi, exval_projforb/4)
[fieldt3, fieldw3, psi3, relE3, conv3, niter3, mallniterc3, Jterms3, maxgrad3, alpha3, invHess3] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb01, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw2, @(w) 1e4*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options2, 13*pi, 2000, 7, 7, 1e-6);
options2 = optionsOCqn(1e-8, 1e4);
options2.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);
options2.invHess0 = invHess2;
[fieldt3, fieldw3, psi3, relE3, conv3, niter3, mallniterc3, Jterms3, maxgrad3, alpha3, invHess3] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, fJforb01, Hoperations, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldw2, @(w) 1e4*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options2, 13*pi, 2000, 7, 7, 1e-8);
Jterms2
Jterms3
figure
plot(0:13*pi/2e3:13*pi, fieldt2)
hold on
plot(0:13*pi/2e3:13*pi, fieldt3)
max(abs(fieldt2-fieldt3))
max(abs(fieldt2-fieldt3), [], 2)
figure
plot(0:1/13:200/13, fieldw3(:,1:201))
save pcontrol2qubits
figure
U3 = reshape(psi3, 125,4);
U3 = reshape(psi3, [125,4]);
size(psi3)
figure
plot(conj(psi3).*psi3)
plot(conj(psi3(1:125,:)).*psi3(1:125,:))
plot(0:13*pi/2e3:13*pi, conj(psi3(1:125,:)).*psi3(1:125,:))
plot(0:13*pi/2e3:13*pi, conj(psi3(126:250,:)).*psi3(126:250,:))
plot(0:13*pi/2e3:13*pi, conj(psi3(251:375,:)).*psi3(251:375,:))
plot(0:13*pi/2e3:13*pi, conj(psi3(376:500,:)).*psi3(376:500,:))
[Jmax_overlap, phases] = Uoverlap_gate(psi3(:, end), Utarget(:), [125 250 375])
plot(0:13*pi/2e3:13*pi, conj(psi3(1:125,:)).*psi3(1:125,:))
cell(2,3)
cl = cell(2,3)
cl{1} = zeros(3,1)
cl{2} = zeros(2,7)
cl{2,3} = zeros(3,7)
whos
clear cl
cl = cell(3)
clear cl
cl = cell(2,3)
cl{2,3} = zeros(3,7)
cl{2,3} = ones(1,7)
clear cl
figure
plot(0:13*pi/2e3:13*pi, fieldt3)
Enl
whos
Enlg
doc speye
cl = cell(2,3)
cl{2,3} = ones(1,7)
cl{2,3}(5)
[sqpsih, proj_vecsh] = projs_parent(psi3, [5 5 5], 1);
[sqpsih00, proj_vecsh] = projs_parent(psi3(1:125,:), [5 5 5], 1);
proj_vecsh
figure
plot(0:13*pi/2e3:13*pi, sqpsih00)
sqpsih01 = projs_parent(psi3(126:250,:), [5 5 5], 1);
figure
plot(0:13*pi/2e3:13*pi, sqpsih00)
plot(0:13*pi/2e3:13*pi, sqpsih01)
max(sqpsih01, [], 2)
sqpsih10 = projs_parent(psi3(251:375,:), [5 5 5], 1);
figure
plot(0:13*pi/2e3:13*pi, sqpsih10)
max(sqpsih10, [], 2)
sqpsih11 = projs_parent(psi3(376:500,:), [5 5 5], 1);
plot(0:13*pi/2e3:13*pi, sqpsih11)
max(sqpsih11, [], 2)
[sqpsiqa00, proj_vecsqa] = projs_parent(psi3(1:125,:), [5 5 5], 1);
max(sqpsiqa00, [], 2)
proj_vecsqa
[sqpsiqa00, proj_vecsqa] = projs_parent(psi3(1:125,:), [5 5 5], 2);
proj_vecsqa
max(sqpsiqa00, [], 2)
sqpsiqa01 = projs_parent(psi3(126:250,:), [5 5 5], 2);
figure
plot(0:13*pi/2e3:13*pi, sqpsiqa01)
max(sqpsiqa01, [], 2)
sqpsiqa10 = projs_parent(psi3(251:375,:), [5 5 5], 2);
max(sqpsiqa10, [], 2)
figure
plot(0:13*pi/2e3:13*pi, sqpsiqa10)
sqpsiqa11 = projs_parent(psi3(375:500,:), [5 5 5], 2);
sqpsiqa11 = projs_parent(psi3(376:500,:), [5 5 5], 2);
figure
plot(0:13*pi/2e3:13*pi, sqpsiqa11)
2*pi/4
ans/2
save pcontrol2qubits
%-- 20/02/2021 21:41 --%
load pcontrol2qubits
sigmax = [0 1; 1 0]
[transferM, Dx] = eig(sigmax)
transferM'
U04x = [kron(transferM(:,1), transferM(:,1)), kron(transferM(:,1), transferM(:,2)), kron(transferM(:,2), transferM(:,1)), kron(transferM(:,2), transferM(:,2))]
transferM4 = kron(transferM, transferM)
transferM4-U04
transferM4-U04x
transferM4^2
transferM^2
iswap = [1 0 0 0; 0 0 1i 0; 0 1i 0 0; 0 0 0 1]
iswapx = transferM4*iswap*transferM4
beep
U0
U0x = sparse(125, 4);
U0x([1 2 6 7], :) = transferM4
transferM4
Utargetx = sparse(125, 4);
Utargetx([1 2 6 7], :) = iswapx
iswapx
Hoperations.psi = @(psi, params, v) H05u*v - params(1)*(Hd1u*v) - params(2)*(Hd2u*v) + params(3)*(Hcplusu*v) + conj(params(3))*(Hcminusu*v);
Hoperationsx.psi = @(psi, params, v) H05u*v - (params(1)+0.4)*(Hd1u*v) - (params(2)+0.4)*(Hd2u*v) + params(3)*(Hcplusu*v) + conj(params(3))*(Hcminusu*v);
Hoperationsx.chi = Hoperationsx.psi
Hoperationsx.diff_chi = Hoperations.diff_chi
Hoperationsx.diff_psi = Hoperations.diff_psi = @(psi1, params1, psi2, params2) (params2(1) - params1(1, :)).*(Hd1u*psi1) + (params2(2) - params1(2, :)).*(Hd2u*psi1)...
+ (params1(3, :) - params2(3)).*(Hcplusu*psi1) + conj(params1(3, :) - params2(3)).*(Hcminusu*psi1);
Hoperations.diff_psi = @(psi1, params1, psi2, params2) (params2(1) - params1(1, :)).*(Hd1u*psi1) + (params2(2) - params1(2, :)).*(Hd2u*psi1)...
+ (params1(3, :) - params2(3)).*(Hcplusu*psi1) + conj(params1(3, :) - params2(3)).*(Hcminusu*psi1);
Hoperationsx.diff_psi = @(psi1, params1, psi2, params2) (params2(1) - params1(1, :)).*(Hd1u*psi1) + (params2(2) - params1(2, :)).*(Hd2u*psi1)...
+ (params1(3, :) - params2(3)).*(Hcplusu*psi1) + conj(params1(3, :) - params2(3)).*(Hcminusu*psi1)
Hoperationsx.diff_chi = Hoperationsx.diff_psi;
[fieldtx, fieldwx, psix, relEx, convx, niterx, mallnitercx, Jtermsx, maxgradx, alphax, invHessx] = OClimf_gate1(U0x(:), Utargetx(:), (1:3)*125, [], Hoperationsx, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e2*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
figure
plot(0:13*pi/2e3:13*pi, fieldtx)
plot(0:13*pi/2e3:13*pi, fieldtx+0.4)
Jtermsx
Jforbx = fJforb_eq(psix, 13*pi, penalforbvu)
[fieldtx1, fieldwx1, psix1, relEx1, convx1, niterx1, mallnitercx1, Jtermsx1, maxgradx1, alphax1, invHessx1] = OClimf_gate1(U0x(:), Utargetx(:), (1:3)*125, [], Hoperationsx, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], fieldwx, @(w) 10*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
Jtermsx1
figure
plot(0:13*pi/2e3:13*pi, fieldtx1+0.4)
Jforbx1 = fJforb_eq(psix1, 13*pi, penalforbvu)
Jforbx1/(13*pi)
figure
plot(0:1/13:200/13, fieldwx(:,1:201))
plot(0:1/13:200/13, fieldwx1(:,1:201))
[fieldtdc, fieldwdc, psidc, relEdc, convdc, niterdc, mallnitercdc, Jtermsdc, maxgraddc, alphadc, invHessdc] = OClimf_gate1(U0(:), Utarget(:), (1:3)*125, [], Hoperationsx, @(t) exp(1i*13*t), 4, fcouplingOp, [-60 10], @(w) 0.3*sin(2*pi*w)*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), @(w) 1e2*0.5.*(1-tanh(2*(w-13))).*heaviside(15-w), options, 13*pi, 2000, 7, 7, 1e-6);
Jtermsdc
Jterms
Jtermsx
figure
plot(0:13*pi/2e3:13*pi, fieldtdc+0.4)