U0 = eye(3)
U0target = [0 1 0; 1 0 0; 0 0 0]
H0
H0u = kron(eye(3), H0)
sigmax = [0 1; 1 0]
whos
miu
miuu = kron(eye(3), miu)
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget, H0, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, T, 0.05, 7, 7, 1e-4, 1e3);
Utarget = [0 1 0; 1 0 0; 0 0 0]
clear U0target
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget, H0, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, T, 0.05, 7, 7, 1e-4, 1e3);
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, T, 0.05, 7, 7, 1e-4, 1e3);
T
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.05, 7, 7, 1e-4, 1e3);
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e4*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.05, 7, 7, 1e-4, 1e3);
Edomain
H0
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.025, 7, 7, 1e-4, 1e3);
mallnitercu
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
mallnitercu
eig(H0u-miuu)
eig(H0u)
Edomain
max(allfieldu)
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
Utarget = [0 1 0; 1 0 0; 0 0 1]
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
options
options.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 30);
[allfieldu, fieldu, psiu, relEu, convu, niteru, mallnitercu, J1u, maxgradu, alphau, invHessu] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
1-J1u
9-J1u
Utarget = [0 1 0; 1 0 0; 0 0 0]
[allfieldu1, fieldu1, psiu1, relEu1, convu1, niteru1, mallnitercu1, J1u1, maxgradu1, alphau1, invHessu1] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, @(t) 0.1*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))).*sin(2*pi*5*t), @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
Utarget = [0 1 0; 1 0 0; 0 0 1]
figure
psiu(:, end).*conj(psiu(:,end))
psiu(:, end)
[allfieldu1, fieldu1, psiu1, relEu1, convu1, niteru1, mallnitercu1, J1u1, maxgradu1, alphau1, invHessu1] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, allfieldu1, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
Utarget = [0 1 0; 1 0 0; 0 0 0]
[allfieldu1, fieldu1, psiu1, relEu1, convu1, niteru1, mallnitercu1, J1u1, maxgradu1, alphau1, invHessu1] = OCqnTDpenal(U0(:), Utarget(:), H0u, Edomain, miuu, allfieldu1, @(t) 1e3*0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))), options, 10, 0.01, 7, 7, 1e-4, 1e3);
1-J1u1
4-J1u1
figure
plot(0:0.1:10, fieldu)
plot(0:0.01:10, fieldu)
hold on
plot(0:0.01:10, fieldu1)
psiu1(:,end).*conj(psiu1(:,end))
psiu(:,end).*conj(psiu(:,end))
figure
plot(0:0.01:10, conj(psiu1(1:3,:).*psiu1(1:3,:))
plot(0:0.01:10, conj(psiu1(1:3,:).*psiu1(1:3,:)))
plot(0:0.01:10, conj(psiu1(1:3,:)).*psiu1(1:3,:))
figure
plot(0:0.01:10, conj(psiu1(4:6,:)).*psiu1(4:6,:))
hold on
plot(0:0.01:10, conj(psiu1(4:6,:)).*psiu1(4:6,:))
clf
plot(0:0.1:10, fieldu1)
plot(0:0.01:10, fieldu1)
hold on
plot(10:-0.01:0, fieldu1)
max(abs(fieldu1-fieldu1(end:-1:1)))
4-J1u1
figure
plot(0:0.01:10, conj(psiu1(1:3,:)).*psiu1(1:3,:))
save unitary3LS allfieldu1 fieldu1 psiu1 relEu1 convu1 niteru1 mallnitercu1 J1u1 maxgradu1 alphau1 H0u targetU U0
whos
save unitary3LS allfieldu1 fieldu1 psiu1 relEu1 convu1 niteru1 mallnitercu1 J1u1 maxgradu1 alphau1 U0 Utarget H0u Edomain miuu