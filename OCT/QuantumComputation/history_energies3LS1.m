0.2/5
2-ans
1/(2*pi*5)
energies3LSn
10/(10*pi)
energies3LSn
conv(niter+1)
figure
plot((0:0.01:3)*5, allfield(1:6:end))
1-J1
Epenal*sum(integw.'.*allfield.^2)
allJ1
energies3LSn
allJ1
whos
allJ1
allconv(:,end)
energies3LS
allconv(:,end)
clear all
energies3LS
energies3LSn
allJ1
1-allJ1
energies
allconv(:,end)
size(allconv)
[allfields(4, 1:(Nti*(Nt_ts - 1) + 1)), fields(4, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), 4), ~, convi, allniter(4), ~, allJ1(4)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.01*sin(t), 1e-3, options, allT(4), dt, 7, 7, 1e-4, 1e4);
allniter
[allfields(4, 1:(Nti*(Nt_ts - 1) + 1)), fields(4, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), 4), ~, convi, allniter(4), ~, allJ1(4)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.01*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
allniter
allJ1
1-allJ1
[allfield7, fields(4, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), 4), ~, convi, allniter(4), ~, allJ1(4)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.01*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
[allfield7, fields(4, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), 4), ~, convi, allniter(4), ~, allJ1(4)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.005*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
[allfield7, fields(4, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), 4), ~, convi, allniter(4), ~, allJ1(4)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.03*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
[allfield7, fields(4, 1:(Nti + 1)), allpsi(:, 1:(Nti + 1), 4), ~, convi, allniter(4), ~, allJ1(4)] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.02*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
[allfield7, field7, allpsi7, ~, conv7, niter7, ~, J17] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.02*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
conv7(end)
J17
1-J17
1e3*(J17-conv7(end))
allconv
1e3*(J17-conv7(end))
energies
allJ1
1-allJ1
energies(4) = 1e3*(J17-conv7(end))
allJ1(4) = J17
1-allJ1
figure
plot((4:10)*5, energies)
max_fields = max(abs(allfields,[],2))
max_fields = max(abs(allfields),[],2)
max_field7 = max(abs(allfield7))
max_fields(4) = max(abs(allfield7))
pi./max_fields
(pi./max_fields)./(2*pi*5*(4:10))
(pi./max_fields)./(2*pi*5*(4:10).')
mean_power = energies./((4:10).'*2*pi*5)
mean_power = energies./((4:10)*2*pi*5)
sq_mean_power = sqrt(mean_power)
clear sq_mean_power
rms_fields = sqrt(mean_power)
pi./(rms_fields*sqrt(2))
pi./(rms_fields*sqrt(2))./allT
T_pi_pulse = pi./(rms_fields*sqrt(2))./allT
fracT = allT./T_pi_pulse
T_pi_pulse = pi./(rms_fields*sqrt(2))
fracT = allT./T_pi_pulse
figure
plot((0:0.01:7)*5, allfield7(1:8:end))
figure
plot((4:10)*5, fracT)
energiesTLSn
save energies3LSn
Api_pulse = pi./allT
energies_pi_pulse = allT.*Api_pulse.^2/2
energies./energies_pi_pulse
figure
fracE = energies./energies_pi_pulse
plot((4:10)*5, fracE)
save energies3LSn
energiesTLSn
energies2./energies_pi_pulse
save energies3LSn
energiesTLSn
energies2./energies_pi_pulse
energiesTLSn
energies2./energies_pi_pulse
1-allJ1_2
energies2/1e2
energiesTLSn
1-allJ1_2
figure
plot((4:10)*5, fracE)
plot((4:10)*5, fracT)
energies_pi_pulse
pi^2./(2*allT)
energies3LSn_longT
allTa
energiesa
pi^2./(2*allTa)
energiesa(1:2)./(pi^2./(2*allTa(1:2)))
energies3LSn_shortT
pi^2./(2*allTb)
energies3LSn_shortT
[allfield7, field7, allpsi7, ~, conv7, niter7, ~, J17] = OCqn(psi0, target, H0, Edomain, miu, @(t) 0.02*sin(t), 1e-3, options, allT(4), dt, 9, 9, 1e-5, 1e4);
energies3LSn_shortT
1-allJ1a
1-allJ1b
energiesb./(pi^2./(2*allTb))
allenergies = [energiesb, energies, energiesa]
allenergies = [energiesb, energies, energiesa(1:2)]
allTtot = [allTb, allT, allTa(1:2)]
allenergies./(pi^2./(2*allTtot))
figure
fracEtot = allenergies./(pi^2./(2*allTtot)
fracEtot = allenergies./(pi^2./(2*allTtot))
allenergies_pi_pulse = (pi^2./(2*allTtot)
allenergies_pi_pulse = (pi^2./(2*allTtot))
allenergies./allenergies_pi_pulse
plot(allTtot, fracEtot)
plot((allTtot)/(5*2*pi), fracEtot)
max(allpsib(3,:,1).*conj(allpsib(3,:,1)))
max(allpsib(3,:,2).*conj(allpsib(3,:,2)))
max(allpsib(3,:,3).*conj(allpsib(3,:,3)))
for j=1:7, max(allpsi(3,:,j).*conj(allpsi(3,:,j))), end
for j=1:2, max(allpsia(3,:,j).*conj(allpsia(3,:,j))), end
plot((allTtot)/(2*pi), fracEtot)
figure
plot(log((allTtot)/(2*pi)), log(fracEtot))
save energies3LSn
xlabel('$f_0t$', 'interpreter', 'latex')
xlabel('$f_0T$', 'interpreter', 'latex')
ylabel('$E/E_\pi$',  'interpreter', 'latex')
plot([0, 150], [1, 1], '-')
hold on
plot([0, 150], [1, 1], '-')
ylabel('$E(t)/E_\pi(T)$',  'interpreter', 'latex')
ylabel('$E(T)/E_\pi(T)$',  'interpreter', 'latex')
plot(log((allTtot)/(2*pi)), log(fracEtot-1))
polyfit(log((allTtot)/(2*pi)), log(fracEtot-1), 1)
xlabel('$\ln(f_0T)$', 'interpreter', 'latex')
ylabel('$\ln(E(T)/E_\pi(T)-1)$',  'interpreter', 'latex')
ylabel('$\ln[E(T)/E_\pi(T)-1]$',  'interpreter', 'latex')
figure
plot(log((allTtot)/(2*pi)), log(allenergies))
polyfit(log((allTtot)/(2*pi)), log(allenergies), 1)
mdl = fitlm(log((allTtot)/(2*pi)),log(fracEtot-1))
figure
whos
3*5/1000
1/ans
Nt_maxb
round(1e3/3)
dtb
dtb/(2*pi)
333*ans
figure
plot(0:dtb/(2*pi):333*dtb/(2*pi), fieldsb(1, 1:334))
plot(0:dtb/(2*pi):667*dtb/(2*pi), fieldsb(2, 1:668))
plot(0:dtb/(2*pi):3/(2*pi), fieldsb(2, 1:668))
plot(0:dtb/(2*pi):3/(2*pi), fieldsb(3, :))
plot(0:dtb/(2*pi):1e3*dtb/(2*pi), fieldsb(3, :))
H0e = diag([0 1 2])
miu*H0e - Hoe*miu
miu*H0e - H0e*miu
miu
M1 = miu*H0e - H0e*miu
miu*M1-M1*miu
M2 = miu*M1-M1*miu
M3 = miu*M2-M2*miu
[allfielde, fielde, allpsie, ~, conve, nitere, ~, J1e] = OCqn(psi0, target, H0e, Edomain, miu, @(t) pi/allT(7)*sin(t), 1e-3, options, allT(7), dt, 7, 7, 1e-5, 1e4);
optionse = optionsOCqn(1e-4, 1e4);
optionse.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 1);
[allfielde, fielde, allpsie, ~, conve, nitere, ~, J1e] = OCqn(psi0, target, H0e, Edomain, miu, @(t) pi/allT(7)*sin(t), 1e-3, optionse, allT(7), dt, 7, 7, 1e-5, 1e4);
optionse.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
[allfielde, fielde, allpsie, ~, conve, nitere, ~, J1e] = OCqn(psi0, target, H0e, Edomain, miu, @(t) pi/allT(7)*sin(t), 1e-3, optionse, allT(7), dt, 7, 7, 1e-5, 1e4);
[allfielde, fielde, allpsie, ~, conve, nitere, ~, J1e] = OCqn(psi0, target, H0e, Edomain, miu, @(t) pi/allT(7)*sin(t), 1e-3, optionse, allT(7), dt/10, 9, 9, 1e-5, 1e4);
[allfielde, fielde, allpsie, ~, conve, nitere, ~, J1e] = OCqn(psi0, target, H0e, Edomain, miu, @(t) pi/allT(7)*sin(t), 1e-3, optionse, allT(7), dt, 9, 9, 1e-5, 1e4);
figure
plot((0:dt:allT(7))/(2*pi), fielde)
J1e
1-J1e
optionse1 = optionsOCqn(1e-4, 1e5);
optionse1.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 10);
[allfielde1, fielde1, allpsie1, relEe1, conve1, nitere1, mallniterce1, J1e1, maxgrade1, alphae1, invHesse1] = OCqn(psi0, target, H0e, Edomain, miu, allfielde, 1e-3, optionse, allT(7), dt, 9, 9, 1e-5, 1e4);
[allfielde1, fielde1, allpsie1, relEe1, conve1, nitere1, mallniterce1, J1e1, maxgrade1, alphae1, invHesse1] = OCqn(psi0, target, H0e, Edomain, miu, allfielde.', 1e-3, optionse, allT(7), dt, 9, 9, 1e-5, 1e4);
pi/allT(7)
[allfielde1, fielde1, allpsie1, relEe1, conve1, nitere1, mallniterce1, J1e1, maxgrade1, alphae1, invHesse1] = OCqn(psi0, target, H0e, Edomain, miu, @(t) 0.2*sin(t), 1e-3, optionse, allT(7), dt, 9, 9, 1e-5, 1e4);
1-J1e
[allfielde2, fielde2, allpsie2, relEe2, conve2, nitere2, mallniterce2, J1e2, maxgrade2, alphae2, invHesse2] = OCqn(psi0, target, H0e, Edomain, miu, allfielde1, 1e-3, optionse1, allT(7), dt, 9, 9, 1e-5, 1e4);
[allfielde2, fielde2, allpsie2, relEe2, conve2, nitere2, mallniterce2, J1e2, maxgrade2, alphae2, invHesse2] = OCqn(psi0, target, H0e, Edomain, miu, allfielde1.', 1e-3, optionse1, allT(7), dt, 9, 9, 1e-5, 1e4);
optionse1.invHess0 = invHesse1;
[allfielde2, fielde2, allpsie2, relEe2, conve2, nitere2, mallniterce2, J1e2, maxgrade2, alphae2, invHesse2] = OCqn(psi0, target, H0e, Edomain, miu, allfielde1.', 1e-3, optionse1, allT(7), dt, 9, 9, 1e-5, 1e4);
J1e2
1-J1e2
figure
plot((0:dt:allT(7))/(2*pi), fielde)
hold on
plot((0:dt:allT(7))/(2*pi), fielde1)
plot((0:dt:allT(7))/(2*pi), fielde2)
figure
plot((0:dt:allT(7))/(2*pi), fielde2)
J1e-conve2(end)
allJ1
energies
(J1e-conve2(end))*1e3
figure
plot((0:dt:allT(7))/(2*pi), conj(psie2).*psie2)
plot((0:dt:allT(7))/(2*pi), conj(allpsie2(1:8:end)).*psie2(1:8:end))
plot((0:dt:allT(7))/(2*pi), conj(allpsie2(1:8:end)).*allpsie2(1:8:end))
plot((0:dt:allT(7))/(2*pi), conj(allpsie2(:,1:8:end)).*allpsie2(:,1:8:end))
size(allpsie2(:,1:8:end))
size(allpsie2(:,1:9:end))
size(allpsie2(:,1:8:end))
size((0:dt:allT(7)))
size(allpsie2)
plot((0:dt:allT(7))/(2*pi), conj(allpsie2).*allpsie2)
whos
clear optionse1 invHesse1 invHesse2
save energies3LSn