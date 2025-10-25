figure
plot(0:1e-4:0.3, field03)
hold on
plot(0.3, energy03, '*')
energy03 = J103-conv03(end)
plot(0.3, energy03, '*')
plot(0.3, energy03*1e4, '*')
figure
plot([0.5, 1:10], log10(1-allJ105), '-o')
hold on
plot(0.3, log10(1-J103), '*', 'r')
plot(0.3, log10(1-J103), '*')
xlabel('$T$ (ns)', 'interpreter', 'latex')
ylabel('$\log_{10}\left(1-\left<\phi|psi(T)\right>\right)$',  'interpreter', 'latex')
ylabel('$\log_{10}\left(1-\left<\phi|\psi(T)\right>\right)$',  'interpreter', 'latex')
psi0E = zeros(7,1);
targetE = zeros(7,1);
psi0E(1) = 1;
targetE(2) = 1;
H0f = diag(E32(1:7));
miuf = XE32(1:7, 1:7);
[allfieldf, fieldf, psif, relEf, convf, niterf, mallnitercf, J1f, maxgradf, alphaf, invHessf] = OCqn(psi0E, targetE, H0f, [0, 166], miuf, @(t) sin(10*pi*t), 1e-3, [], 4, 1e-2, 7, 7, 1e-4, 1e3);
mallnitercf
allJ1f
J1f
allJ1
max(allconv, [], 2)
figure
plot(0:0.01:4, conj(psif).*psif)
[allfieldf3, fieldf3, psif3, relEf3, convf3, niterf3, mallnitercf3, J1f3, maxgradf3, alphaf3, invHessf3] = OCqn(psi0E(1:3), targetE(1:3), H0f(1:3, 1:3), [0, , miuf(1:3, 1:3), @(t) sin(10*pi*t), 1e-3, [], 4, 1e-2, 7, 7, 1e-4, 1e3);
E(1:3)
[allfieldf3, fieldf3, psif3, relEf3, convf3, niterf3, mallnitercf3, J1f3, maxgradf3, alphaf3, invHessf3] = OCqn(psi0E(1:3), targetE(1:3), H0f(1:3, 1:3), [0, 81] , miuf(1:3, 1:3), @(t) sin(10*pi*t), 1e-3, [], 4, 1e-2, 7, 7, 1e-4, 1e3);
J1f3
figure
plot(0:0.01:4, conj(psif3).*psif3)
figure
plot(0:0.01:4, fieldf)
hold on
plot(0:0.01:4, fieldf3)
[allfieldf1, fieldf1, psif1, relEf1, convf1, niterf1, mallnitercf1, J1f1, maxgradf1, alphaf1, invHessf1] = OCqn(psi0E, targetE, H0f, [0, 166], miuf, allfieldf3, 1e-3, [], 4, 1e-2, 7, 7, 1e-4, 1e3);
[allfieldf1, fieldf1, psif1, relEf1, convf1, niterf1, mallnitercf1, J1f1, maxgradf1, alphaf1, invHessf1] = OCqn(psi0E, targetE, H0f, [0, 166], miuf, allfieldf3.', 1e-3, [], 4, 1e-2, 7, 7, 1e-4, 1e3);
J1f1
J1f
J1f3
plot(0:0.01:4, fieldf1)
figure
plot(0:0.01:4, conj(psif1).*psif1)
figure
plot(t, 0.5*(tanh(2*(t-2.5))-tanh(2*(t-7.5))))
xlabel('$t$ (ns)', 'interpreter', 'latex')
ylabel('$s(t)$',  'interpreter', 'latex')
figure
plot(allt_grid, allfield)
xlabel('$t$ (ns)', 'interpreter', 'latex')
yl
ylabel('$\epsilon(t)$',  'interpreter', 'latex')
J1
1-J1
conv(end)
H0
H0/(2*pi)
figure
whos
energies05
plot([0.3 0.5, 1:10], energies03)
plot([0.5, 1:10], energies05)
plot(1./[0.5, 1:10], energies05)
plot(log10(1./[0.5, 1:10]), log10(energies05))
fieldw = dctIfrom_ig_sym1(allfield, 10, 0.05*(1-cos((0:5)*pi/6))/2, w(1:6)*0.05);
figure
plot(0:pi/10:pi/0.05, fieldw)
max(abs(fieldt-dctI(fieldw)))
max(abs(field-dctI(fieldw)))
figure
plot(0:0.05:10, dctI(fieldw))
hold on
plot(0:0.05:10, field)
plot(0:0.05:10, field./dctI(fieldw))
fieldw = dctIfrom_ig_sym1(allfield, 10, 0.05*(1-cos((0:5)*pi/6))/2, w(1:6));
max(abs(field-dctI(fieldw)))
max(abs(field-dctI(fieldw)./field))
max(abs(allfield-dctIint_grid(fieldw)./field))
max(abs(allfield-dctIintgrid(fieldw)./field))
max(abs(allfield-dctIintgrid(fieldw,10, )./field))
max(abs(allfield-dctIintgrid(fieldw,10, 0.05*(1-cos((0:5)*pi/6))/2)./field))
max(abs(allfield-dctIintgrid(fieldw,10, 0.05*(1-cos((0:5)*pi/6))/2)))
max(abs(allfield-dctIintgrid(fieldw,10, 0.05*(1-cos((0:5)*pi/6))/2)./allfield))
max(abs(allfield-dctIintgrid(fieldw,10, 0.05*(1-cos((0:5)*pi/6))/2)))
figure
plot(allt_grid, allfield)
hold on
plot(allt_grid, dctIintgrid(fieldw,10, 0.05*(1-cos((0:5)*pi/6))/2)))
plot(allt_grid, dctIintgrid(fieldw,10, 0.05*(1-cos((0:5)*pi/6))/2))
size(fieldw)
noise = rand(1,201)-0.5;
max(noise)
min(noise)
noise = (rand(1,201)-0.5)*2;
min(noise)
max(noise)
field_noise = dctIintgrid(fieldw.*(1+noise*1e-2),10, 0.05*(1-cos((0:5)*pi/6))/2);
figure
plot(allt_grid, allfield)
hold on
plot(allt_grid, field_noise)
max(allfield-field_noise)
max(abs(allfield-field_noise))
max(abs((allfield-field_noise)./allfield))
figure
plot(all_tgrid, allfield-field_noise)
plot(allt_grid, allfield-field_noise)
max(abs(allfield-field_noise))/max(abs(allfield))
field_noise_tp = dctIintgrid(fieldw.*(1+noise*1e-2), 10, t_ts);
t_ts
dt
tcheb
t_ts_tp = [t_ts, (t_ts(4)+t_ts(5))/2]
t_ts_tp = [t_ts; (t_ts(4)+t_ts(5))/2]
field_noise_tp = dctIintgrid(fieldw.*(1+noise*1e-2), 10, t_ts_tp);
size(field_noise_tp)
t_ts_tp = [t_ts(1:6); (t_ts(4)+t_ts(5))/2]
field_noise_tp = dctIintgrid(fieldw.*(1+noise*1e-2), 10, t_ts_tp);
size(field_noise_tp)
[psi_noise, mniterc, ~, max_errors_psi] = SemiGlobalHparams(@(psi, field, v) H0*v - field*miu*v, @(psi1, field1, psi2, field2) -miu*(psi1*diag(field1 - field2)), 1, field_noise_tp, [], [], psi0, 0:0.05:10, 200, 7, 7, 1e-8, 10, 16, [], false);
[psi_noise, mniterc, ~, max_errors_psi] = SemiGlobalHparams(@(psi, field, v) H0*v - field*miu*v, @(psi1, field1, psi2, field2) -miu*(psi1*diag(field1 - field2)), 1, field_noise_tp, [], [], psi0, 0:0.05:10, 2000, 7, 7, 1e-8, 10, 16, [], false);
[psi_noise, mniterc, ~, max_errors_psi] = SemiGlobalHparams(@(psi, field, v) H0*v - field*miu*v, @(psi1, field1, psi2, field2) -miu*(psi1*diag(field1 - field2)), 1, field_noise_tp, [], [], psi0, 0:0.05:10, 200, 7, 7, 1e-8, 10, 16, [], false);
[psi_noise, mniterc, ~, max_errors_psi] = SemiGlobalHparams(@(psi, field, v) H0*v - field*miu*v, @(psi1, field1, psi2, field2) -miu*(psi1*diag(field1 - field2)), 1, field_noise_tp, [], EdomainGHz, psi0, 0:0.05:10, 200, 7, 7, 1e-8, 10, 16, [], false);
psi_noise(end).*conj(psi_noise(end))
size(psi_noise)
psi_noise(:,end).*conj(psi_noise(:,end))
1-psi_noise(:,end).*conj(psi_noise(:,end))
noise_effect
fidelity
noise_effect
fidelity
figure
noise_effect
fidelity
field_noise1 = dctIintgrid(fieldw.*(1+noise*1e-1), 10, t_ts(1:6));
noise_effect
fidelity
field_noisei[7:7:end] = [];
field_noisei(7:7:end) = [];
field_noise_tpi(7:7:end) = [];
size(field_noise_tpi)
size(field_noise1)
max(abs(field_noise1-field_noise_tpi))
figure
plot(allt_grid, allfield)
hold on
plot(allt_grid, field_noise_tpi)
figure
plot(allt_grid, field_noise_tpi-allfield)
mean(abs(field_noise_tpi-allfield))
ans./max(abs(field_noise_tpi-allfield))
max(abs(field_noise_tpi-allfield))
ans./max(abs(field_noise_tpi-allfield))
mean(abs(field_noise_tpi-allfield))
ans./max(abs(field_noise))
fidelity
(3+6+4+4+6+4+5+4+6+5   +(45+30+30+40+30+30+20)/60)*8.5
75.5*8.5
(3+6+4+4+6+4+5+4+6+5+4+3+6+5+5   +(45+30+30+40+30+30+20+30+40+20+15)/60)*8.5
A = [1 2 3; 4 5 6; 7 8 9]
v = [1 2 3]
A.*v
v.*A
A.*v.'
C = zeros(2, 2, 2)
v(1:2).'*C
v(1:2).'.*C
v(1:2).*C
C(:,:,1) = A(1:2,1:2)
C(:,:,2) = A(2:3,2:3)
v(1:2).'.*C
v(1:2).*C
A(1:2,1:2)*C