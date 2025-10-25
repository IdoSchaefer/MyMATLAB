sigmax = [0 1; 1 0]
H0 = [0 0; 0 1];
H04 = kron(eye(2), H0)
mu4 = kron(eye(2), sigmax)
U0 = eye(2)
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(U0(:), sigmax(:), H04, [-1 3], mu4, @(t) sin(t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqnTDpenal(U0(:), sigmax(:), H04, [-1 3], mu4, @(t) sin(t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
psi(:,end).*conj(psi(:,end))
J1
1-J1
4-J1
mallniterc
figure
plot(0:0.1:10, field)
[psi01, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, [-1 2], sigmax, [1; 0], [0 T], 100, 7, 7, 1e-5, allfield);
[psi01, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, [-1 2], sigmax, [1; 0], [0 10], 100, 7, 7, 1e-5, allfield);
psi01(:,end).*conj(psi01(:,end))
[psi10, ~, mniterc] = solveOCMLSih1(@ihfieldMLS, H0, [-1 2], sigmax, [0; 1], [0 10], 100, 7, 7, 1e-5, allfield);
psi10(:,end).*conj(psi10(:,end))
figure
plot(0:0.1:10, psi01.*conj(psi01)
plot(0:0.1:10, psi01.*conj(psi01))
plot(0:0.1:10, psi01(1:6:end).*conj(psi01(1:6:end)))
size(psi01)
plot(0:0.1:10, psi01(1:6:end,:).*conj(psi01(1:6:end,:)))
plot(0:0.1:10, psi01(:,1:6:end).*conj(psi01(:,1:6:end)))
hold on
plot(0:0.1:10, psi10(:,1:6:end).*conj(psi10(:,1:6:end)))
[allfieldy, fieldy, U, relEy, convy, nitery, mallnitercy, J1y, maxgrady, alphay, invHessy] = OCqnTDpenal(U0(:), sigmay(:), H04, [-1 2], mu4, @(t) sin(t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
sigmay = [0 1i; -1i 0]
sigmax*sigmay
sigmay = [0 -1i; 1i 0]
sigmax*sigmay
[allfieldy, fieldy, U, relEy, convy, nitery, mallnitercy, J1y, maxgrady, alphay, invHessy] = OCqnTDpenal(U0(:), sigmay(:), H04, [-1 2], mu4, @(t) sin(t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
[allfieldy, fieldy, U, relEy, convy, nitery, mallnitercy, J1y, maxgrady, alphay, invHessy] = OCqnTDpenal(U0(:), sigmay(:), H04, [-1 2], mu4, @(t) sin(t).*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
[allfieldy, fieldy, U, relEy, convy, nitery, mallnitercy, J1y, maxgrady, alphay, invHessy] = OCqnTDpenal(U0(:), sigmay(:), H04, [-1 2], mu4, @(t) sin(t)*0.5.*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 30);
[allfieldy, fieldy, U, relEy, convy, nitery, mallnitercy, J1y, maxgrady, alphay, invHessy] = OCqnTDpenal(U0(:), sigmay(:), H04, [-1 2], mu4, @(t) sin(t)*0.5.*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
J1y
1-J1y
4-J1y
figure
plot(0:0.1:10, fieldy)
plot(0:0.1:10, conj(U).*U)
U0(:)
U(:,end)
[0 1 1 0]*U(:,end)
[0 1].'.*U(:,1:2)
[0 1]*U(:,1:2)
[0 1]*U(1:2,end)
[1 0]*U(3:4,end)
[0 1]*psi(1:2,end)
[1 0]*psi(3:4,end)
[allfield1, field1, U1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqnTDpenal(U0(:), [0;1;0;1], H04, [-1 2], mu4, @(t) 0.1*sin(t)*0.5.*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.1, 7, 7, 1e-4, 1e3);
J11
1-J11
2-J11
U(:,end).*conj(U(:,end))
whos
clear allfield1 field1 U1 relE1 conv1 niter1 mallniterc1 J11 maxgrad1 alpha1 invHess1 invHess invHessy
save unitaryTLS