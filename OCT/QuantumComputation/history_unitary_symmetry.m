load('unitary3LS.mat')
miu = [0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0]
miuy = [0 -1i 0; 1i 0 -sqrt(2)*1i; 0 sqrt(2)*1i 0]
miuyu = kron(eye(3), miuy);
[allfielduy, fielduy, psiuy, relEuy, convuy, niteruy, mallnitercuy, J1uy, maxgraduy, alphauy, invHessuy] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuyu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
options = optionsOCqn(1e-4, 1e4);
[allfielduy, fielduy, psiuy, relEuy, convuy, niteruy, mallnitercuy, J1uy, maxgraduy, alphauy, invHessuy] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuyu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 30);
[allfielduy, fielduy, psiuy, relEuy, convuy, niteruy, mallnitercuy, J1uy, maxgraduy, alphauy, invHessuy] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuyu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
[allfielduy, fielduy, psiuy, relEuy, convuy, niteruy, mallnitercuy, J1uy, maxgraduy, alphauy, invHessuy] = OCqnTDpenal(U0(1:6), Utarget(1:6), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
U0
U0(1:6)
U0(:)
U0((1:6).')
[allfielduy, fielduy, psiuy, relEuy, convuy, niteruy, mallnitercuy, J1uy, maxgraduy, alphauy, invHessuy] = OCqnTDpenal(U0((1:6).'), Utarget((1:6).'), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
1-J1uy/4
figure
plot(0:0.01:10, fielduy)
[allfielduy1, fielduy1, psiuy1, relEuy1, convuy1, niteruy1, mallnitercuy1, J1uy1, maxgraduy1, alphauy1, invHessuy1] = OCqnTDpenal(U0((1:6).'), Utarget((1:6).'), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), allfieldy, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
[allfielduy1, fielduy1, psiuy1, relEuy1, convuy1, niteruy1, mallnitercuy1, J1uy1, maxgraduy1, alphauy1, invHessuy1] = OCqnTDpenal(U0((1:6).'), Utarget((1:6).'), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), allfielduy, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e4);
figure
plot(0:0.01:10, fielduy1)
1-J1uy1/4
options1 = optionsOCqn(1e-5, 1e4);
options1.invHess0 = invHessuy;
[allfielduy1, fielduy1, psiuy1, relEuy1, convuy1, niteruy1, mallnitercuy1, J1uy1, maxgraduy1, alphauy1, invHessuy1] = OCqnTDpenal(U0((1:6).'), Utarget((1:6).'), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), allfielduy, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options1, 10, 0.01, 7, 7, 1e-5, 1e4);
1-J1uy1/4
figure
plot(0:0.01:10, fielduy1)
options2 = optionsOCqn(1e-8, 1e4);
options2.invHess0 = invHessuy1;
[allfielduy2, fielduy2, psiuy2, relEuy2, convuy2, niteruy2, mallnitercuy2, J1uy2, maxgraduy2, alphauy2, invHessuy2] = OCqnTDpenal(U0((1:6).'), Utarget((1:6).'), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), allfielduy1, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options2, 10, 0.01, 7, 7, 1e-5, 1e4);
figure
plot(0:0.01:10, fielduy2)
max(abs(fielduy2+fielduy2(end:-1:1)))
Utarget
Utargety = [0 -1i 0; 1i 0 0]
Utargety = [0 -1i; 1i 0; 0 0]
[allfielduty, fielduty, psiuty, relEuty, convuty, niteruty, mallnitercuty, J1uty, maxgraduty, alphauty, invHessuty] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options2, 10, 0.01, 7, 7, 1e-8, 1e4);
1-J1uty/4
figure
plot(0:0.01:10, fielduty)
[allfielduyty, fielduyty, psiuyty, relEuyty, convuyty, niteruyty, mallnitercuyty, J1uyty, maxgraduyty, alphauyty, invHessuyty] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options2, 10, 0.01, 7, 7, 1e-8, 1e4);
options3 = optionsOCqn(1e-8, 1e4);
options3.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 30);
[allfielduty, fielduty, psiuty, relEuty, convuty, niteruty, mallnitercuty, J1uty, maxgraduty, alphauty, invHessuty] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options3, 10, 0.01, 7, 7, 1e-8, 1e4);
1-J1uty/4
figure
plot(0:0.01:10, fielduty)
[allfielduyty, fielduyty, psiuyty, relEuyty, convuyty, niteruyty, mallnitercuyty, J1uyty, maxgraduyty, alphauyty, invHessuyty] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options3, 10, 0.01, 7, 7, 1e-8, 1e4);
options4 = optionsOCqn(1e-8, 1e4);
options4.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 40);
[allfielduyty, fielduyty, psiuyty, relEuyty, convuyty, niteruyty, mallnitercuyty, J1uyty, maxgraduyty, alphauyty, invHessuyty] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options4, 10, 0.01, 7, 7, 1e-8, 1e4);
options4.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 100);
[allfielduyty, fielduyty, psiuyty, relEuyty, convuyty, niteruyty, mallnitercuyty, J1uyty, maxgraduyty, alphauyty, invHessuyty] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options4, 10, 0.01, 7, 7, 1e-8, 1e4);
options4 = options3;
[allfielduyty1, fielduyty1, psiuyty1, relEuyty1, convuyty1, niteruyty1, mallnitercuyty1, J1uyty1, maxgraduyty1, alphauyty1, invHessuyty1] = OCqnTDpenal(U0((1:6).'), Utargety(:), H0u(1:6, 1:6), Edomain, miuyu(1:6,1:6), allfielduyty, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options3, 10, 0.01, 7, 7, 1e-8, 1e4);
1-J1uyty/4
J1uyty
1-J1uyty1/4
figure
plot(0:0.01:10, fielduyty1)
figure
plot(0:0.01:10, fieldu1)
hold on
plot(0:0.01:10, fielduy2)
figure
plot(0:0.01:10, fielduty)
hold on
plot(0:0.01:10, fielduyty1)
whos
save unitary_symmetry allfieldu1 allfielduty allfielduy2 allfielduyty1 fieldu1 fielduty fielduy2 fielduyty1