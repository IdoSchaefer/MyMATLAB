Vf = @(x) 1 - 1./sqrt(x.^2+1)
40/64
[Vabs480, xabs480, x480, p480, K480, Nx480, V0480] = get_prop_vars(Vf, [-480 480], 40/64, Vopt649x);
figure
plot(x480, real(Vabs480))
plot(x480, imag(Vabs480))
save coulomb_optV480 Vabs480 xabs480 x480 p480 K480 Nx480 V0480 Vf Vopt649x
load('fields_13_15_17.mat')
[fieldtgn, fieldwgn, psign, evatgn, evawgn, evmiutgn, evmiuwgn, mnitercgn, Jgn, J1gn, J2gn, Jorthgn, Jpnormgn] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015)*sum(fieldt.^2)/sum(fieldtg.^2), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4
[fi0480, E0480, ~, E480, ~, ~] = gsV(Vf, [-480 480],
480/0.625
960/0.625
[fi0480, E0480, ~, E480, ~, ~] = gsV(Vf, [-480 480], 1536);
save coulomb_optV480 Vabs480 xabs480 x480 p480 K480 Nx480 V0480 Vf Vopt649x fi0480 E0480 E480
[fieldt1b480, fieldw1b480, psi1b480, evat1b480, evaw1b480, evmiut1b480, evmiuw1b480, mniterc1b480, J1b480, J11b480, J21b480, Jorth1b480, Jpnorm1b480] = guessresults_pnaE0b(fi0480, Vabs480, 1, [-480 480], xabs480, -x480./(1 + x480.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw1b, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
plot(w(1:1001), evaw1b(1:1001)))
plot(w(1:1001), evaw1b(1:1001))
hold on
plot(w(1:1001), evaw1b480(1:1001))
figure
plot(w(1:1001), evaw1b(1:1001) - evaw1b480(1:1001))
plot(w(1:1001), evaw1b(1:1001) - evaw1b480(1:1001).')
plot(w(1:1001), (evaw1b(1:1001) - evaw1b480(1:1001).')./evaw1b(1:1001))
figure
plot(w, log10(abs(evaw1b)))
hold on
plot(w, log10(abs(evaw1b480)))
plot(w(1:1001), evaw1b(1:1001) - evaw1b480(1:1001).')
[fieldt7a480, fieldw7a480, psi7a480, evat7a480, evaw7a480, evmiut7a480, evmiuw7a480, mniterc7a480, J7a480, J17a480, J27a480, Jorth7a480, Jpnorm7a480] = guessresults_pnaE0b(fi0480, Vabs480, 1, [-480 480], xabs480, -x480./(1 + x480.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw7a, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.84).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
plot(w, log10(abs(evaw7a)))
hold on
plot(w, log10(abs(evaw7a480)))
figure
plot(w(1:1001), evaw1b(1:1001))
hold on
clf
plot(w(1:1001), evaw7a(1:1001))
hold on
plot(w(1:1001), evaw7a480(1:1001))
[fieldt3a480, fieldw3a480, psi3a480, evat3a480, evaw3a480, evmiut3a480, evmiuw3a480, mniterc3a480, J3a480, J13a480, J23a480, Jorth3a480, Jpnorm3a480] = guessresults_pnaE0b(fi0480, Vabs480, 1, [-480 480], xabs480, -x480./(1 + x480.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw3a, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
max(abs(eva1b-eva1b480))
max(abs(evaw1b-evaw1b480))
max(abs(evaw1b-evaw1b480.'))
max(abs(evaw1b-evaw1b480.'))/1.81
max(abs(evaw7a-evaw7a480.'))/1.692
figure
plot(w(1:1001), evaw3a(1:1001))
hold on
plot(w(1:1001), evaw3a480(1:1001))
figure
plot(w, log10(abs(evaw3a)))
hold on
plot(w, log10(abs(evaw3a480)))
max(abs(evaw3a-evaw3a480.'))/0.9028
figure
plot(w(1:1001), evaw3a - eva3a480)
plot(w(1:1001), evaw3a - evaw3a480)
plot(w, evaw3a - evaw3a480)
plot(w, evaw3a - evaw3a480.')
plot(w(1:1001), evaw3a(1:1001) - evaw3a480(1:1001).')
[fieldt6480, fieldw6480, psi6480, evat6480, evaw6480, evmiut6480, evmiuw6480, mniterc6480, J6480, J16480, J26480, Jorth6480, Jpnorm6480] = guessresults_pnaE0b(fi0480, Vabs480, 1, [-480 480], xabs480, -x480./(1 + x480.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), fieldw6, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-1.02).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
plot(w(1:1001), evaw6(1:1001))
hold on
plot(w(1:1001), evaw6480(1:1001))
max(abs(evaw6-evaw6480.'))/0.9894
figure
plot(w(1:1001), evaw6(1:1001) - evaw6480(1:1001).')
figure
plot(w(1:1001), evaw1b(1:1001) - evaw1b480(1:1001).')
plot(w(1:1001), evaw7a(1:1001) - evaw7a480(1:1001).')
plot(w(1:1001), evaw3a(1:1001) - evaw3a480(1:1001).')
(J11b - J11b480)
(J11b - J11b480)/J11b
(J13a - J13a480)/J13a
(J17a - J17a480)/J17a
(J16 - J16480)/J16
save evaw480 evaw1b480 evaw3a480 evaw7a480 evaw6480 J11b480 J13a480 J17a480 J16480