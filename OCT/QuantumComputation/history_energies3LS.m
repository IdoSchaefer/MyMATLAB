energies3LS
energies
figure
plot(1:10, energies)
[allfield6, field6 psi6, relE6, conv6, niter6, mallniterc6, J16, maxgrad6, alpha6, invHess6] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, options, 6, 0.01, 7, 7, 1e-4, 1e3);
[allfield6, field6 psi6, relE6, conv6, niter6, mallniterc6, J16, maxgrad6, alpha6, invHess6] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(5*t), 1e-3, options, 6, 0.01, 7, 7, 1e-4, 1e3);
allJ1
allJ1(6) = J16
1-allJ1
energies(6) = 1e3*(allJ1(6) - conv6(end))
[allfield1, field1 psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(5*t), 1e-3, options, 1, 0.005, 7, 7, 1e-4, 1e3);
options1 = optionsOCqn(1e-4, 1e3);
options1.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 20);
[allfield1, field1 psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(5*t), 1e-3, options1, 1, 0.005, 7, 7, 1e-4, 1e3);
options1.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 30);
options1.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 40);
[allfield1, field1 psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(5*t), 1e-3, options1, 1, 0.0025, 7, 7, 1e-4, 1e3);
[allfield1, field1 psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, options1, 1, 0.0025, 7, 7, 1e-4, 1e3);
allJ11-conv1(end)
J11-conv1(end)
energies
(J11-conv1(end))*1e3
[allfield1, field1 psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(pi*t), 1e-3, options1, 1, 0.0025, 7, 7, 1e-4, 1e3);
(J11-conv1(end))*1e3
max(allfields, 2)
max(allfields, 1)
max(allfields)
doc max
max(allfields, [], 2)
max(abs(allfields, [], 2))
max(abs(allfields), [], 2))
max(abs(allfields), [], 2)
max(abs(allfield1))
figure
plot(0:1e-3:1, field1)
size(field1)
plot(0:0.0025:1, field1)
J11
allJ1(1) = J11
1-allJ1
[allfield05, field05, psi05, relE05, conv05, niter05, mallniterc05, J105, maxgrad05, alpha05, invHess05] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(pi*t), 1e-3, options1, 0.5, 0.001, 7, 7, 1e-4, 1e3);
J105
(J105-conv05(end))*1e3
[allfield05, field05, psi05, relE05, conv05, niter05, mallniterc05, J105, maxgrad05, alpha05, invHess05] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, options1, 0.5, 0.001, 7, 7, 1e-4, 1e3);
J105
(J105-conv05(end))*1e3
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, options1, 0.025, 0.001, 7, 7, 1e-4, 1e3);
options1 = optionsOCqn(1e-4, 1e3);
options1.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 40);
options2 = optionsOCqn(1e-4, 1e3);
options2.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 80);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, options2, 0.025, 0.001, 7, 7, 1e-4, 1e3);
J1025
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, options2, 0.025, 0.001, 7, 7, 1e-4, 1e3);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, options1, 0.025, 0.001, 7, 7, 1e-4, 1e3);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, options2, 0.025, 0.001, 7, 7, 1e-4, 1e3);
options2.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 60);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, options2, 0.025, 0.001, 7, 7, 1e-4, 1e3);
options2.f_max_alpha =  @(field, direction) alpha_max_x(field, direction, 100);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, options2, 0.025, 0.001, 7, 7, 1e-4, 1e3);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, [], 0.025, 0.001, 7, 7, 1e-4, 1e3);
J1025
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-5, [], 0.025, 0.001, 7, 7, 1e-4, 1e3);
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-5, [], 0.025, 0.0005, 7, 7, 1e-4, 1e3);
J105
J1025
[allfield025, field025, psi025, relE025, conv025, niter025, mallniterc025, J1025, maxgrad025, alpha025, invHess025] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.025, 0.0005, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, [], 0.03, 0.0005, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, [], 0.03, 0.0005, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-5, [], 0.03, 0.0005, 7, 7, 1e-4, 1e3);
J103
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 0.0005, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 0.0001, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 0.00001, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 1e-5, 7, 7, 1e-4, 1e3);
doc round
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 1e-5, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 1e-6, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 1e-5, 9, 9, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 5e-6, 9, 9, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-6, [], 0.03, 1e-5, 7, 7, 1e-4, 1e3);
figure
plot(0:0.001:0.5, field05)
energies
energies(1) = 1e3*(J11 - conv1(end))
allJ1
energies05 = [1e3*(J105 - conv05(end)), energies]
allJ105 = [J105, allJ1]
figure
plot([0.5, 1:10], energies05)
figure
plot
plot(0:0.001:0.5, field05)
hold on
plot(0:0.001:1, field1)
plot(0:0.0025:1, field1)
plot(0:0.0025:2, fields(2,:))
plot(0:0.01:2, fields(2,:))
plot(0:0.01:2, fields(2,1:200))
plot(0:0.01:2, fields(2,1:201))
plot(0:0.01:3, fields(3,1:301))
figure
plot(0:0.001:0.5, conj(psi05).*psi05)
figure
plot(0:0.0025:1, conj(psi1).*psi1)
figure
plot(0:0.01:10, conj(allpsi(:,:,10)).*allpsi(:,:,10))
fieldw1 = dctIintgrid(allfield1, 1, 0.0025*(1-cos((0:6)*pi/7))/2);
size(fieldw1)
size(allfield1)
w=chebweights(7,1)
fieldw1 = dctIfrom_ig_sym1(allfield1, 1, 0.0025*(1-cos((0:6)*pi/7))/2, w(1:6));
w(1) = 2*w(1)
fieldw1 = dctIfrom_ig_sym1(allfield1, 1, 0.0025*(1-cos((0:6)*pi/7))/2, w(1:6)*0.0025);
fieldw1 = dctIfrom_ig_sym1(allfield1, 1, 0.0025*(1-cos((0:5)*pi/6))/2, w(1:6)*0.0025);
figure
size(fieldw1)
plot(0:pi:pi/0.0025, fieldw1)
whos
clear invHess03 invHess025 invHess05 invHess1 invHess10 invHess100 invHess6
whos
clear invHess
whos
save energies3LS
%-- 11/06/2020 10:27 --%
whos
(4.6*2)^2
4.6^2
4.6^2/0.2
4.6^2/0.8
26.45/0.4
4.6/0.2
4.6^2/(0.2^2*8)
1.6^2/(0.4*6.626)
8*0.2*2
1/(3.2*pi)
ans*(2*pi*25)^2
(10*pi)^2/(3.2*pi)
ans/pi
(2*4.6*pi)^2/(3.2*pi)
ans/pi
ans/2
ans/4.6
4.6/5*96.6
4.6/0.2
ans/2
2*pi/ans
m=1/(3.2*pi)
[fi0, E0, x, E, P, H] = gsV(@cos, [-pi pi], 12, m);
E
[fi0, E0, x, E, P, H] = gsV(@(x) 1-cos(x), [-pi pi], 12, m);
E
(E(2:12) - E(1:11))/(2*pi)
[fi0, E0, x, E, P, H] = gsV(@(x) 26.45*pi*(1-cos(x)), [-pi pi], 12, m);
E
(E(2:12) - E(1:11))/(2*pi)
figure
plot(x, Vf(x))
Vf = (x) 26.45*pi*(1-cos(x))
Vf = @(x) 26.45*pi*(1-cos(x))
plot(x, Vf(x))
hold on
plot(x, E*ones(1,12))
plot(E, '*')
plot(zeros(12,1), E, '*')
[fi016, E016, x16, E16, P16, H16] = gsV(@(x) 26.45*pi*(1-cos(x)), [-pi pi], 16, m);
E16
plot(x16, Vf(16))
clf
plot(x16, Vf(16))
clf
plot(x16, Vf(x16))
hold on
plot(zeros(12,1), E, '*')
plot(zeros(16,1), E16, 'o')
E16(1:12) - E
E16
(E16(1:12) - E)./(E16(1:12)
(E16(1:12) - E)./(E16(1:12))
(E16(2:7) - E16(1:6))/(2*pi*4.6)
E16(1)/(2*pi*4.6)
E16(1:7)/(2*pi*4.6)
E16(2:7)-E16(1:6)
(E16(2:7)-E16(1:6))/(2*pi)
[fi032, E032, x32, E32, P32, H32] = gsV(@(x) 26.45*pi*(1-cos(x)), [-pi pi], 32, m);
E32(1:8)
E32(1:7)-E16(1:7)
(E32(1:7)-E16(1:7))./E32(1:7)
(E16(2)-E16(1))
(E16(2)-E16(1))/(2*pi)
(E16(3) - 2*E16(2)+E16(1))/(2*pi)
20*ans
XE16 = P16\diag(x16)*P16;
XE16
XE16(1:7, 1:7)
XE32 = P32\diag(x32)*P32;
XE32(1:7, 1:7)
E32(1:7).'
figure
plot([0.5, 1:10], energies05)
hold on
xlabel('$T$ (ns)', 'interpreter', 'latex')
ylabel('$\int_0^T\epsilon^2(t)\,dt$')
ylabel('$\int_0^T\epsilon^2(t)\,dt$',  'interpreter', 'latex')
[allfieldf, fieldf, psif, relEf, convf, niterf, mallnitercf, J1f, maxgradf, alphaf, invHessf] = OCqn(psi0f, targetf, H0f, Edomainf, miuf, @(t) sin(10*pi*t), 1e-3, [], 4, 1e-2, 7, 7, 1e-4, 1e3);
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-3, [], 0.3, 0.0001, 7, 7, 1e-4, 1e3);
energies03 = [1e3*(J103 - conv03(end)), energies05]
J103
mallniterc03
[allfield03, field03, psi03, relE03, conv03, niter03, mallniterc03, J103, maxgrad03, alpha03, invHess03] = OCqn(psi0, target, H0, Edomain, miu, @(t) sin(10*pi*t), 1e-4, [], 0.3, 0.0001, 7, 7, 1e-4, 1e3);
J103
J103-conv03(end)
figure
clear invHess03
whos