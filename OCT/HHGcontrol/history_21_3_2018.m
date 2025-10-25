whos
Vf
figure
plot(x, Vf(x))
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, m, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
options = optionsOCqn(1e-3, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
(minf - f0)/(ro*-df0_dalpha)
(minf - f0)/(ro*df0_dalpha)
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
figure
plot(0:0.05:100, allfield(1:8:allt_lasti))
options.max_x = 5;
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
options.max_x = 10;
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
a0M
figure
plot(0:0.05:100, allfield(1:8:allt_lasti))
allapsi(:,1)
chiihterm(:,1)
mallniterc
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
mallniterc
alpha
figure
J1
w=0:pi/100:pi/0.05;
plot(w(1:101), fieldw(1:101))
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
J1
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
mallniterc
maxgrad
figure
plot(0:0.05:100, fieldt)
plot(w(1:101), fieldw(1:101))
figure
plot(w(1:1001), evmiuw(1:1001))
plot(w(1:101), evmiuw(1:101))
plot(w(1:121), evmiuw(1:121))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
figure
plot(0:niter, conv)
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, [fieldw, zeros(1,2000)], @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.025, 9, 9, 1e-3, 1e3
options.invHess0 = invHess;
options.alpha0 = alpha
options.alpha0 = 3.8551e-03
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, [fieldw, zeros(1,2000)], @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.025, 9, 9, 1e-3, 1e3);
options = optionsOCqn(1e-3, 1e3);
options.max_x = 10;
maxgrad1
maxgrad
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, [fieldw, zeros(1,2000)], @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.025, 9, 9, 1e-3, 1e3);
options.alpha0 = 3.8551e-03;
options.invHess0 = invHess;
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, [fieldw, zeros(1,2000)], @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.025, 9, 9, 1e-3, 1e3);
eigvalH=diag(eig(invHess))
eigvalH=(eig(invHess))
options = optionsOCqn(1e-3, 1e3);
options.max_x = 10;
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[falpha, Dfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/10:b);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/10:b);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/100:b);
hold on
allDalphaFD = (allfalpha(2:end) - allfalpha(1:end-1))/(Lsection/100);
plot(a:Lsection/100:Lsection*99/100, allDfalphaFD)
plot(a:Lsection/100:Lsection*99/100, allDalphaFD)
Lsection/100
allDfalpha(1:end-1)./allDalphaFD
doc diff
allDfalphaFD = gradient(allfalpha, Lsection/100);
plot(a:Lsection/100:Lsection*99/100, allDfalphaFD)
plot(a:Lsection/100:b, allDfalphaFD)
allDfalpha./allDfalphaFD
plot(a:Lsection/100:b, allDfalphaFD*30)
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/100:b);
allDfalphaFD = gradient(allfalpha, Lsection/100);
allDfalpha./allDfalphaFD
hold on
plot(a:Lsection/100:b, allDfalphaFD*30)
[allfalphag, allDfalphag] = sample_falpha(fDf, x0, , a:Lsection/100:b
[allfalphag, allDfalphag] = sample_falpha(fDf, x0, -fgrad0, a:Lsection/100:b);
1/dw
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], miu, @(t) 1, 0.1, [], 10, 0.05, 5, 5, 1e-3, 1e3
options = optionsOCqn(1e-3, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], miu, @(t) 1, 0.1, [], 10, 0.05, 5, 5, 1e-3, 1e3
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], miu, @(t) 1, 0.1, [], 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, [], 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
options.alpha0 = 0.5;
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
options.alpha0 = 0.9;
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
options.alpha0 = 1e-2;
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
options.alpha0 = 1;
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, [], 10, 0.05, 5, 5, 1e-3, 1e3);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, 0:0.01:1);
allDfalphaFD = gradient(allfalpha, 0.01);
hold on
plot(0:0.01:1, allDfalphaFD)
options.alpha0 = 2;
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
maxgrad
options.alpha0 = 1;
options.invHess0 = invHess;
options = optionsOCqn(1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], allfield, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], allfield.', 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], allfield.'+(rand(801,1)-0.5)*0.1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], allfield.'+(rand(801,1)-0.5)*0.5, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], allfield.'+(rand(801,1)-0.5)*2, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 2, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
miu = [0 1 1;
1 0 1;
1 1 0];
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, 0:0.01:1);
allDfalphaFD = gradient(allfalpha, 0.01);
hold on
plot(0:0.01:1, allDfalphaFD)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2).*cos(10*w), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, 0:0.01:1);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/100:b);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, b:Lsection/100:a);
allDfalphaFD = gradient(allfalpha, 0.01);
hold on
plot(0:0.01:1, allDfalphaFD)
dw
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 10*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
options.max_x = 10;
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/100:b);
allDfalphaFD = gradient(allfalpha, 0.01);
hold on
plot(a:Lsection/100:b, allDfalphaFD)
allDfalphaFD = gradient(allfalpha, Lsection/100);
plot(a:Lsection/100:b, allDfalphaFD)
figure
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 1000*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
figure
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 1000*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
figure
plot(0:dt:T, dctI(fieldw))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 1000*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2).*cos(10*w), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
figure
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 1000*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
figure
plot(w(1:121), evmiuw(1:121))
options.alpha_estimation = true
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 1000*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt1, fieldw1, psi1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) rectanglefun(w, 0, 1.3), @(w) 1000*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
figure
plot(0:niter, conv)
hold on
plot(0:niter1, conv1)
figure
plot(w(1:121), evmiuw(1:121))
hold on
plot(w(1:121), evmiuw1(1:121))
J1
J11
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2).*cos(10*w), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
[fieldt2, fieldw2, psi2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, fieldw1, @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3
options.invHess0 = invHess1;
options.invHess0 = [];
options1=options;
options1.invHess = invHess1;
options1.alpha0 = alpha1;
[fieldt2, fieldw2, psi2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, fieldw1, @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options1, 100, 0.05, 9, 9, 1e-3, 1e3);
[fieldt2, fieldw2, psi2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, fieldw1, @(w) rectanglefun(w, 0, 1.3), @(w) 100*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
J12
J11
conv1(end)-J11
conv1(end)/10-J11
(conv1(end)-J11)+J11/10
J12-conv2(end)
conv2(end)-J12
(conv2(end)-J12)+J12/10
(conv2(end)-J12)+J12/2
J12/(conv2(end)-J12)
[fieldt3, fieldw3, psi3, evmiut3, evmiuw3, relE3, conv3, niter3, mallniterc3, J13, maxgrad3, alpha3, invHess3] = OCfxE0b_miu_qn(fi0, Vf, 1, [-10 85], xdomain, x, fieldw2, @(w) rectanglefun(w, 0, 1.3), @(w) 100/1.5*rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
J13
J12
figure
plot(w(1:121), evmiuw2(1:121))
plot(w(1:121), evmiuw3(1:121))
plot(w(1:121), fieldw3(1:121))
plot(0:0.05:100, fieldt)
plot(0:0.05:100, fieldt3)
plot(0:0.05:100, fieldt2)
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
options.alpha_estimation = false;
[allfield, field, psi, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
78/31
86/38
load aviv1
whos
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
options = optionsOCqn(1e-3, 1e3);
options.alpha_estimation = true;
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 7, 7, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 5e2, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3
options.max_x = 2*max(abs(fieldg))
fieldg
options.max_x = 0.002;
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
options.alpha_estimation = false;
options.alpha_estimation = true;
options.max_x = 0.01;
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
options.max_x = 0.005;
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
whos
figure
plot(0:niter, conv)
plot(0:niter1, conv1)
Epenal
options.alpha_estimation = false;
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg*0.1, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
fieldtg
fieldg
fieldg1 = @(t)1e-4.*exp(-0.00000000005*(t-5e5).^2).*cos(-0.000000001*(t-5e5).^2)
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg1, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg1, Epenal, [], T, 1e3, 5, 5, 1e-3, 1e3);
options.alpha_estimation = true;
[allfield1, field1, psi1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCqn_op(psi0all, projO_all, H0all, [-1e-4, 1e-3], miu_all, fieldg1, Epenal, options, T, 1e3, 5, 5, 1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
figure
plot(0:0.05:100, dctI([minusgrad; zeros(2001-42, 0)])/dctfactor)
plot(0:0.05:100, dctI([minusgrad; zeros(2001-42, 1)])/dctfactor)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
options.f_termination(dif_f, dif_x, fmin_old, gradmin_old, xmin_old)
options = optionsOCqn(1e-3, 1e3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
figure
plot(0:0.05:100, dctI([minusgrad; zeros(2001-42, 1)])/dctfactor)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
allfield(end)
plot(0:0.05:100, dctI([minusgrad; zeros(2001-42, 1)])/dctfactor)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
chebweights(5, dt)/dt
itegw2 = ans
ans*2
igweights
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
sqrt(2/pi)*sum([0.5*vfilterE(1); vfilterE(2:Nt); 0.5*vfilterE(Nt + 1)])*dw
ans - dctfilterE0
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
dctfilterET - sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1); vfilterE(2:Nt).*coswT(2:Nt); 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
max(abs(J1fun - 0.5*evmiuw.^2.*vfiltermiu.*integw))
max(abs(J2fun +(fieldwnz.^2./vfilterEnz).*integwnz))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
max(abs(minusgrad - 2*integwnz.*(fieldwnz./vfilterEnz - dct_chimiupsi(iEnz) + lambda0 + lambdaT*coswT(iEnz))))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
fieldt_unc0 - sqrt(2/pi)*sum(fieldw_unc(iEnz).*integw)
fieldt_unc0 - sqrt(2/pi)*sum(fieldw_unc(iEnz).'.*integw)
fieldt_unc0 - sqrt(2/pi)*sum(fieldw_unc(iEnz).*integwnz)
fieldt_unc0 - sqrt(2/pi)*sum(fieldw_unc(iEnz).*integw(iEnz))
fieldt_unc0
fieldt_uncT
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
fieldt_unc0
fieldt_uncT
[allfield(1), allfield( nnnnnnnnend)]
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
figure
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
dctfilterE0 - sqrt(2/pi)*sum(vfilterEnz.*integwnz)
dctfilterET - sqrt(2/pi)*sum(vfilterEnz.*coswT(iEnz).*integwnz)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfE0b_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) rectanglefun(w, 0, 1.3).*cos(10*w), @(w) 0.25*rectanglefun(w, 0, 1.3), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
figure
plot(0:0.05:100, fieldt)
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
Vf
whos
xdomain
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0, Vf, 1, xdomain, x, exp(-x) - 1, @(x) 0, @(x) 0, @(w) rectanglefun(w, 0, 1.3).*cos(w*5*2*pi/1.3), @(w) 100*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 3.3, 3.5), options, 100, 0.05, 9, 9, 1e-3, 1e3);
niter/nprop
niter/412
412/niter
figure
plot(0:0.05:100, fieldt)
w=0:pi/100:pi/0.05;
plot(w(1:101), fieldw(1:101))
options
plot(w(1:1001), evaw(1:1001))
plot(w(1:1001), evmiuw(1:1001))
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options.max_x = 15;
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options.max_x = 10;
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
penalnormf(normpsiT)
sum(J1fun) + sum(J2fun)
sum(J1fun)
sum(J2fun)
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/100:b);
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, b:Lsection/5:a);
allDfalphaFD = gradient(allfalpha, Lsection/5);
hold on, plot(b:Lsection/5:a, allDfalphaFD)
fprintf('\nWarning: The alpha1 value (%d) was limited by the user preferences (iteration No. %d).\n', 1, 2)
figure
plot(w(1:1001), eva(1:1001))
plot(w(1:1001), evaw(1:1001))
w=0:pi/1000:pi/0.2;
plot(w(1:1001), evaw(1:1001))
psi(:, end).*conj(psi(:,end))
sum(psi(:, end).*conj(psi(:,end)))
figure
plot(0:0.2:1e3, fieldt)
dctfactor = 1e3/(sqrt(5e3*pi))
hold on
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015)))
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
plot(0:0.2:1e3, fieldt)
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 10*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldt1, fieldw1, psi1, evat1, evaw1, evmiut1, evmiuw1, relE1, conv1, niter1, mallniterc1, J11, maxgrad1, alpha1, invHess1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 10*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
sum(psi(:, end).*conj(psi(:,end)))
sum(psi1(:, end).*conj(psi1(:,end)))
figure
plot(w(1:101), fieldw(1:101))
hold on
plot(w(1:101), fieldw1(1:101))
plot(w(1:101), 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.015))
J11
figure
plot(w(1:1001), evaw(1:1001))
plot(w(1:1001), log10(evaw(1:1001)))
figure
plot(0:0.2:1e3, fieldt1)
hold on
plot(0:0.2:1e3, fieldt)
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
2*pi/0.06
iw = instwcos(fieldt, 1e3);
figure
plot(0:0.2:1e3, iw)
iw1 = instwcos(fieldt1, 1e3);
plot(0:0.2:1e3, iw1)
options.max_x = 15;
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 10*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options.invHess0 = invHess1;
options.alpha0 = alpha1;
[fieldt2, fieldw2, psi2, evat2, evaw2, evmiut2, evmiuw2, relE2, conv2, niter2, mallniterc2, J12, maxgrad2, alpha2, invHess2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw1, @(w) 10*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
184/78
J12
figure
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
hold on
plot(0:0.2:1e3, fieldt2)
figure
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.015))
hold on
plot(w(1:101), fieldw2(1:101))
iw2 = instwcos(fieldt2, 1e3);
figure
plot(0:0.2:1e3, iw2)
figure
plot(w(1:1001),
plot(w(1:1001), log10(evaw2(1:1001)))
figure
plot(w(1:1001), (evaw2(1:1001)))
maxgrad2
maxgrad1
maxgrad
sum(psi1(:, end).*conj(psi1(:,end)))
sum(psi2(:, end).*conj(psi2(:,end)))
doc colormap
Wfield = psi2wigner(fieldt(1:10:end));
Wfield = psi2wigner(fieldt(1:10:end).');
mesh(0:1e3, 0:5*pi/1e3:pi/0.2, Wfield)
size(Wfield)
size(fieldt(1:10:end))
Wfield2 = psi2wigner(fieldt2(1:10:end-1).');
size(Wfield)
clear Wfield
size(Wfield2)
mesh(0:999, 0:5*pi/1e3:(pi/0.2-5*pi/1e3), Wfield2)
mesh(0:999, -pi/10:5*pi/1e3:(pi/10-5*pi/1e3), Wfield2(:, 400:599))
mesh(0:999, -pi/10:5*pi/1e3:(pi/10-5*pi/1e3), Wfield2(:, 300:699))
size(-pi/10:5*pi/1e3:(pi/10-5*pi/1e3))
mesh(0:999, -pi/10:pi/1e3:(pi/10-pi/1e3), Wfield2(:, 400:599))
size(-pi/10:pi/1e3:(pi/10-pi/1e3))
size(500:599)
size(400:599)
size(Wfield2)
mesh(-pi/10:pi/1e3:(pi/10-pi/1e3), 0:999, Wfield2(:, 400:599))
mesh(0:999, -pi/10:pi/1e3:(pi/10-pi/1e3), Wfield2(:, 400:599).')
mesh(0:999, -pi/10:pi/1e3:(pi/10-pi/1e3), Wfield2(:, 400:599))
surf(0:999, -pi/10:pi/1e3:(pi/10-pi/1e3), Wfield2(:, 400:599), 'LineStyle', 'none')
mesh(0:999, -pi/10:pi/1e3:(pi/10-pi/1e3), Wfield2(400:599, :))
Weva2 = psi2wigner(evat2(1:10:end-1).');
figure
mesh(0:999, -pi/10:pi/1e3:(pi/10-pi/1e3), Weva2(400:599, :))
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), Weva2(400:599, :))
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), Weva2)
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), log10(Weva2))
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), real(log10(Weva2)))
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), Weva2)
viewP(Weva2, -pi/2:pi/1e3:(pi/2-pi/1e3), 0.1)
viewP(Wfield2, -pi/2:pi/1e3:(pi/2-pi/1e3), 0.1)
viewP(Wfield2, -pi/2:pi/1e3:(pi/2-pi/1e3), 0.01)
502.6/0.2
fieldt2(2514)
fieldt2(2513:2515)
hold on
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)))
figure
plot(0:0.2:1e3, evat2)
328/0.2
evaw2a = dctI(evat2(1:1641));
evaw2b = dctI(evat2(1641:5001));
evaw2a = dctI(evat2(1:1641))*328/(sqrt(1640*pi));
evaw2b = dctI(evat2(1641:5001))*(1e3-328)/sqrt(pi*(5001-1641));
figure
1e3-328
plot(0:pi/328:pi, evaw2a(1:329))
figure
plot(w(1:1001), (evaw2(1:1001)))
figure
plot(0:pi/672:pi, evaw2b(1:672))
plot(0:pi/672:pi, evaw2b(1:673))
figure
plot(w(1:1001), (evaw2(1:1001))/1e3)
hold on
plot(w(1:1001), (evaw2a(1:1001))/328)
plot(w(1:1001), (evaw2b(1:1001))/672)
plot(0:pi/328:pi, evaw2a(1:329)/328)
plot(0:pi/672:pi, evaw2b(1:673)/672)
iweva2 = instwcos(evat2, 1e3);
figure
plot(0:0.2:1e3, iweva2)
0.78/0.57
0.57*13
0.57*15
0.55*15
0.55*14
iweva2 = instw(evat2, 1e3);
figure
plot(0:0.2:1e3, iweva2)
iweva2 = instwcos(evat2, 1e3);
hold on
iweva2b = instw(evat2, 0.2);
plot(0:0.2:1e3, iweva2)
plot(0:0.2:1e3, iweva2b)
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options, 100, 0.05, 5, 5, 1e-3, 1e3);
options1 = optionsOCqn(1e-3, 100);
options1.f_max_alpha = get_f_max_alphaOCf(5, 0.05, 100, @(w) 0.25*exp(-20*(w-1).^2))
options1.f_max_alpha(zeros(54,1), ones(54,1))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options1, 100, 0.05, 5, 5, 1e-3, 1e3);
figure
plot(0:0.05:100, fieldt)
options1.f_max_alpha = get_f_max_alphaOCf(1, 0.05, 100, @(w) 0.25*exp(-20*(w-1).^2))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options1, 100, 0.05, 5, 5, 1e-3, 1e3);
max(allfield)
find(allfield<=1)
find(allfield>1)
find(allfield(1:4:end)>1)
options1.f_max_alpha = get_f_max_alphaOCf(0.5, 0.05, 100, @(w) 0.25*exp(-20*(w-1).^2))
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options1, 100, 0.05, 5, 5, 1e-3, 1e3);
max(allfield)
max(abs(allfield))
figure
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), Weva2)
viewP(Weva2, -pi/2:pi/1e3:(pi/2-pi/1e3), 0.01)
figure
mesh(0:999, -pi/2:pi/1e3:-pi/10, Weva2(1:51,:))
mesh(0:999, -pi/2:pi/1e3:-pi/10, Weva2(1:501,:))
size(-pi/2:pi/1e3:-pi/10)
mesh(0:999, -pi/2:pi/1e3:-pi/10, Weva2(1:401,:))
mesh(0:999, -pi/2:pi/1e3:-pi/5, Weva2(1:201,:))
pi/5
pi/4
3*pi/8
pi/8
mesh(0:999, -pi:2*pi/1e3:pi, Weva2)
mesh(0:999, -pi:2*pi/1e3:(pi-2*pi/1e3), Weva2)
mesh(0:999, -pi/2:pi/1e3:(pi/2-pi/1e3), Weva2)
mesh(0:999, (-pi/2+200*pi/1e3):pi/1e3:(-pi/2+300*pi/1e3), Weva2(201:301,:))
mesh(0:999, (-pi/2+240*pi/1e3):pi/1e3:(-pi/2+260*pi/1e3), Weva2(241:261,:))
figure
mesh(0:999, 0:pi/1e3:pi/10, Wfield2(501:601,:))
mesh(0:999, 10*pi/1e3:pi/1e3:35*pi/1e3, Wfield2(511:536,:))
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options1 = optionsOCqn(1e-3, 1e3);
options1 = optionsOCqn(1e-3, 100);
options = optionsOCqn(1e-3, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
[fieldt, fieldw, psi, evat, evaw, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, alpha, invHess] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
a
b
2*(b - a)*Dfa/(abs(fa) + abs(fb))
(b-a)
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/10:b);
allDfalphaFD = gradient(allfalpha, Lsection/10);
hold on
plot(a:Lsection/10:b, allDfalphaFD)
[allfalpha, allDfalpha] = sample_falpha(fDf, x0, direction, a:Lsection/100:b);
allDfalphaFD = gradient(allfalpha, Lsection/100);
hold on
plot(a:Lsection/10:b, allDfalphaFD)
plot(a:Lsection/100:b, allDfalphaFD)
eig(invHess)
invHess
maxgrad
[fieldta, fieldwa, psia, evata, evawa, evmiuta, evmiuwa, relEa, conva, nitera, mallniterca, J1a, maxgrada, alphaa, invHessa] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldw, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
eig(invHess)
maxgrada
figure
plot(0:0.2:1e3, fieldta)
hold on
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
figure
plot(w(1:101), fieldwa)
plot(w(1:101), fieldwa(1:101))
[fieldtb, fieldwb, psib, evatb, evawb, evmiutb, evmiuwb, relEb, convb, niterb, mallnitercb, J1b, maxgradb, alphab, invHessb] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwa, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
mallnitercb
mallniterca
mallniterc
maxgradb
options1=options
options1.invHess0 = invHessa;
options1 = optionsOCqn(1e-4, 1e3);
options1.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
options1.invHess0 = invHessa;
[fieldtb, fieldwb, psib, evatb, evawb, evmiutb, evmiuwb, relEb, convb, niterb, mallnitercb, J1b, maxgradb, alphab, invHessb] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwa, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options1, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtb, fieldwb, psib, evatb, evawb, evmiutb, evmiuwb, relEb, convb, niterb, mallnitercb, J1b, maxgradb, alphab, invHessb] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwa, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options1, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(w(1:101), fieldwb(1:101))
figure
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
hold on
plot(0:0.3:1e3, fieldtb)
plot(0:0.2:1e3, fieldtb)
figure
plot(w(1:1001), evawb(1:1001))
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwb, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options2 = optionsOCqn(1e-5, 1e3);
options2.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
options2.invHess0 = invHessb;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwb, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options2, 1e3, 0.2, 7, 7, 1e-5, 1e3);
relE
relEb
relEc
maxgradc
maxgradb
relEa
hold on
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)))
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.015))
options2 = optionsOCqn(1e-3, 1e3);
options2.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
options2.maxNiter = 1;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options2, 1e3, 0.2, 7, 7, 1e-3, 1);
[fieldtd, fieldwd, psid, evatd, evawd, evmiutd, evmiuwd, relEd, convd, niterd, mallnitercd, J1d, maxgradd, alphad, invHessd] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options2.maxNiter = 2;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options2, 1e3, 0.2, 7, 7, 1e-3, 2);
[fieldtd, fieldwd, psid, evatd, evawd, evmiutd, evmiuwd, relEd, convd, niterd, mallnitercd, J1d, maxgradd, alphad, invHessd] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options2.maxNiter = 1;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options2, 1e3, 0.2, 7, 7, 1e-3, 1
options2
options2.invHess0 = invHessc;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options2, 1e3, 0.2, 7, 7, 1e-3, 1);
options2.invHess0 = [];
options2.maxNiter = 3;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options2, 1e3, 0.2, 7, 7, 1e-3, 3);
doc fminunc
edit fminunc
edit fminusub
edit lineSearch
figure
plot(0:niterc, convc)
[fieldtd, fieldwd, psid, evatd, evawd, evmiutd, evmiuwd, relEd, convd, niterd, mallnitercd, J1d, maxgradd, alphad, invHessd] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options2.maxNiter = 4;
options3 = optionsOCqn(1e-3, 1e3);
options3.ro = 0.01;
[allfields, fields, psis, relEs, convs, niters, mallnitercs, J1s, maxgrads, alphas, invHesss] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
options3 = optionsOCqn(1e-3, 1e3);
[allfields, fields, psis, relEs, convs, niters, mallnitercs, J1s, maxgrads, alphas, invHesss] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options, 10, 0.05, 5, 5, 1e-3, 1e3);
[allfields, fields, psis, relEs, convs, niters, mallnitercs, J1s, maxgrads, alphas, invHesss] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options3, 10, 0.05, 5, 5, 1e-3, 1e3);
options3.ro = 0.01;
[allfields, fields, psis, relEs, convs, niters, mallnitercs, J1s, maxgrads, alphas, invHesss] = OCqn([1;0], [0;1], diag([1 2]), [-1 4], [0 1; 1 0], @(t) 0.1, 0.1, options3, 10, 0.05, 5, 5, 1e-3, 1e3);
[fieldtt, fieldwt, psit, evmiutt, evmiuwt, relEt, convt, nitert, mallniterct, J1t, maxgradt, alphat, invHesst] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), options3, 100, 0.05, 5, 5, 1e-3, 1e3);
options3 = optionsOCqn(1e-3, 1e3);
options3.ro = 0.01;
[fieldtt, fieldwt, psit, evmiutt, evmiuwt, relEt, convt, nitert, mallniterct, J1t, maxgradt, alphat, invHesst] = OCf_qn([1;0;0], diag([1 1.9 3]), [-1 5], miu, @(w) exp(-20*(w-1).^2), @(w) 0.25*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), [], 100, 0.05, 5, 5, 1e-3, 1e3);
options1
options
options3 = options;
options3.ro = 0.01;
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options3, 1e3, 0.2, 7, 7, 1e-3, 3);
max(abs(fieldtc-fieldt))
J13
J12
J1a
J1b
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 0, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtg, fieldwg, psig, evatg, evawg, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jpnormg] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(w) 0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e7*rectanglefun(w, 0.77, 0.79), 0, 1e3, 0.2, 7, 7, 1e-3);
[Jg, J1g, J2g, Jorthg, Jpnormg]
[fieldtc, fieldwc, psic, evatc, evawc, evmiutc, evmiuwc, relEc, convc, niterc, mallnitercc, J1c, maxgradc, alphac, invHessc] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e7*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldtd, fieldwd, psid, evatd, evawd, evmiutd, evmiuwd, relEd, convd, niterd, mallnitercd, J1d, maxgradd, alphad, invHessd] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3
figure
plot(0:0.2:1e3, fieldtc)
hold on
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
figure
plot(w(1:101), 0.55*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.015))
plot(w(1:101), 0.5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.015))
hold on
plot(w(1:101), fieldwc(1:101))
[fieldtd, fieldwd, psid, evatd, evawd, evmiutd, evmiuwd, relEd, convd, niterd, mallnitercd, J1d, maxgradd, alphad, invHessd] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwc, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
[fieldte, fieldwe, psie, evate, evawe, evmiute, evmiuwe, relEe, conve, nitere, mallniterce, J1e, maxgrade, alphae, invHesse] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwd, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015))/dctfactor)
hold on
plot(0:0.2:1e3, fieldte)
plot(0:0.2:1e3, fieldtb)
plot(0:0.2:1e3, fieldta)
options3 = options;
options3.invHess0 = invHesse;
options3 = optionsOCqn(1e-4, 1e3);
options3.invHess0 = invHesse;
options3.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
[fieldtf, fieldwf, psif, evatf, evawf, evmiutf, evmiuwf, relEf, convf, niterf, mallnitercf, J1f, maxgradf, alphaf, invHessf] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwe, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options3, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(0:0.2:1e3, fieldte)
hold on
plot(0:0.2:1e3, fieldtf)
figure
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.1))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*cos((w-0.06)))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*cos((w-0.06)*pi/0.1))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(5*(w-0.06)*pi))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(20*(w-0.06)*pi))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(1e3*(w-0.06)))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(5e2*(w-0.06)))/dctfactor)
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(7e2*(w-0.06)))/dctfactor)
[fieldth, fieldwh, psih, evath, evawh, evmiuth, evmiuwh, relEh, convh, niterh, mallniterch, J1h, maxgradh, alphah, invHessh] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(7e2*(w-0.06)), @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:0.2:1e3, dctI(0.5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(7e2*(w-0.06)))/dctfactor)
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(7e2*(w-0.06)))/dctfactor)
hold on
plot(0:0.2:1e3, fieldth)
[fieldth1, fieldwh1, psih1, evath1, evawh1, evmiuth1, evmiuwh1, relEh1, convh1, niterh1, mallniterch1, J1h1, maxgradh1, alphah1, invHessh1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:0.2:1e3, fieldth1)
[fieldth2, fieldwh2, psih2, evath2, evawh2, evmiuth2, evmiuwh2, relEh2, convh2, niterh2, mallniterch2, J1h2, maxgradh2, alphah2, invHessh2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh1, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
hold on
plot(0:0.2:1e3, fieldth2)
[fieldth3, fieldwh3, psih3, evath3, evawh3, evmiuth3, evmiuwh3, relEh3, convh3, niterh3, mallniterch3, J1h3, maxgradh3, alphah3, invHessh3] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh2, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:0.2:1e3, fieldth2)
hold on
plot(0:0.2:1e3, fieldth3)
max(abs(fieldth3-fieldth3))
max(abs(fieldth3-fieldth2))
options4 = optionsOCqn(1e-4, 1e3);
options4.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)))
[fieldth3, fieldwh3, psih3, evath3, evawh3, evmiuth3, evmiuwh3, relEh3, convh3, niterh3, mallniterch3, J1h3, maxgradh3, alphah3, invHessh3] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwh2, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options4, 1e3, 0.2, 7, 7, 1e-4, 1e3);
plot(0:0.2:1e3, fieldth3)
figure
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*pi/0.015))
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*700))
hold on
plot(w(1:101), fieldwh3(1:101))
iwh3 = instwcos(fieldth3, 1e3);
figure
plot(0:0.2:1e3, iwh3)
hold on
plot(0:0.2:1e3, fieldth3)
hold on
plot(0:0.2:1e3, fieldth3)
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*700))
plot(0:0.2:1e3, dctI(5*exp(-(w-0.06).^2/(2*0.01^2)).*sin(7e2*(w-0.06)))/dctfactor)
461/0.2+1
fieldth3(2306)
figure
plot(w(1:1001), evawh3(1:1001))
sum(psih3(:, end).*conj(psih3(:,end)))
evawh3a = dctI(fieldth3(1:2306))/sqrt(2305*pi);
evawh3b = dctI(fieldth3(2306:5001))/sqrt((5000 - 2306)*pi);
figure
length(0:pi/461:pi/10)
length(0:pi/539:pi/10)
length(0:pi/539:pi)
length(0:pi/461:pi)
plot(0:pi/461:pi, evawh3(1:462))
evawh3a = dctI(evah3(1:2306))/sqrt(2305*pi);
evawh3a = dctI(evath3(1:2306))/sqrt(2305*pi);
evawh3b = dctI(evath3(2306:5001))/sqrt((5000 - 2306)*pi);
plot(0:pi/461:pi, evawh3(1:462))
plot(0:pi/461:pi, evawh3a(1:462))
hold on
plot(0:pi/461:pi, evawh3b(1:462))
plot(0:pi/539:pi, evawh3b(1:539))
plot(0:pi/539:pi, evawh3b(1:540))
fieldwh3b = dctI(evath3(2306:5001))/sqrt((5000 - 2306)*pi);
fieldwh3b = dctI(fieldth3(2306:5001))/sqrt((5000 - 2306)*pi);
fieldwh3a = dctI(fieldth3(1:2306))/sqrt(2305*pi);
figure
plot(0:pi/539:pi, fieldwh3b(1:540))
plot(0:pi/539:pi/10, fieldwh3b(1:54))
hold on
plot(0:pi/461:pi/10, fieldwh3a(1:46))
plot(0:pi/461:pi/10, fieldwh3a(1:47))
hold on
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*700))
plot(w(1:101), 5*exp(-(w(1:101)-0.06).^2/(2*0.01^2)).*sin((w(1:101)-0.06)*700)/1e3)
78/15
78/11
78/9
78/17
78/19
78/3.4
noise = exp(-(w-0.06).^2/(2*0.01^2)).*(rand(5001)-1);
noise = exp(-(w.'-0.06).^2/(2*0.01^2)).*(rand(5001, 1)-1);
size(noise)
figure
plot(w(1:101), noise(1:101))
noise = exp(-(w.'-0.06).^2/(2*0.01^2)).*(rand(5001, 1)-0.5);
plot(w(1:101), noise(1:101))
figure
plot(0:0.2:1e3, dctI(noise)/dctfactor)
noise_con = fieldw20b(noise, exp(-(w.'-0.06).^2/(2*0.01^2)), pi/1e3);
hold on
plot(w(1:101), noise_con(1:101))
hold on
plot(0:0.2:1e3, dctI(noise_con)/dctfactor)
J12
J1a
J1b
[fieldtb1, fieldwb1, psib1, evatb1, evawb1, evmiutb1, evmiuwb1, relEb1, convb1, niterb1, mallnitercb1, J1b1, maxgradb1, alphab1, invHessb1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwb+noise_con, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
size(fieldwb)
[fieldtb1, fieldwb1, psib1, evatb1, evawb1, evmiutb1, evmiuwb1, relEb1, convb1, niterb1, mallnitercb1, J1b1, maxgradb1, alphab1, invHessb1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwb+noise_con.', @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:niterb, convb)
options5 = options;
options5 = optionsOCqn(1e-4, 1e3);
options5.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)));
options5.invHess0 = invHessb1;
[fieldtb2, fieldwb2, psib2, evatb2, evawb2, evmiutb2, evmiuwb2, relEb2, convb2, niterb2, mallnitercb2, J1b2, maxgradb2, alphab2, invHessb2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwb1, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options5, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(0:0.2:1e3, fieldtb)
hold on
plot(0:0.2:1e3, fieldtb2)
plot(0:0.2:1e3, fieldtb1)
max(abs(fieldtb-fieldtb2))
2*pi/82
2*pi/77
0.06/pi*1e3
pi/1e3
w(20)
cww = zeros(1, 5001);
cww(20) = 20;
figure
plot(0:0.2:1e3, dctI(cww)/dctfactor)
cww_con = fieldw20b(cww, exp(-(w.'-0.06).^2/(2*0.01^2)), pi/1e3);
cww_con = fieldw20b(cww, exp(-(w-0.06).^2/(2*0.01^2)), pi/1e3);
size(w)
size(cww)
cww_con = fieldw20b(cww.', exp(-(w.'-0.06).^2/(2*0.01^2)), pi/1e3).';
hold on
plot(0:0.2:1e3, dctI(cww_con)/dctfactor)
figure
plot(w(1:101), cww_con(1:101))
[fieldti, fieldwi, psii, evati, evawi, evmiuti, evmiuwi, relEi, convi, niteri, mallniterci, J1i, maxgradi, alphai, invHessi] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, cww_con, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
figure
plot(0:0.2:1e3, fieldti)
hold on
plot(0:0.2:1e3, dctI(cww_con)/dctfactor)
figure
plot(w(1:101), cww_con(1:101))
hold on
plot(w(1:101), fieldwi)
plot(w(1:101), fieldwi(1:101))
[fieldti1, fieldwi1, psii1, evati1, evawi1, evmiuti1, evmiuwi1, relEi1, convi1, niteri1, mallniterci1, J1i1, maxgradi1, alphai1, invHessi1] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwi, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options, 1e3, 0.2, 7, 7, 1e-3, 1e3);
options5 = optionsOCqn(1e-4, 1e3);
options5.invHess0 = invHessi1;
options5.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldti2, fieldwi2, psii2, evati2, evawi2, evmiuti2, evmiuwi2, relEi2, convi2, niteri2, mallniterci2, J1i2, maxgradi2, alphai2, invHessi2] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.5*50*sech(50*(x-0.9))^2, fieldwi1, @(w) 1e2*exp(-(w-0.06).^2/(2*0.01^2)), @(w) 4e3*rectanglefun(w, 0.77, 0.79), options5, 1e3, 0.2, 7, 7, 1e-4, 1e3);
figure
plot(0:0.2:1e3, dctI(cww_con)/dctfactor)
hold on
plot(0:0.2:1e3, fieldti)
plot(0:0.2:1e3, fieldti1)
plot(0:0.2:1e3, fieldti2)
max(abs(fieldti1-fieldti2))
max(abs(fieldti1-fieldti2)./abs(fieldti2))
(norm(fieldti1-fieldti2)/norm(fieldti2))
figure
plot(w(1:101), cww_con(1:101))
hold on
plot(w(1:101), fieldwi(1:101))
plot(w(1:101), fieldwi2(1:101))
figure
plot(w(1:1001), evawi2(1:1001))
iwi2 = instwcos(fieldti2, 1e3);
figure
plot(0:0.2:1e3, iwi2)
figure
plot(0:0.2:1e3, fieldti2)
iwi2 = instwcos([fieldti2, fieldti2(end-1:-1:1)], 1e3);
plot(0:0.2:2e3, iwi2)
iwi2 = instwcos([fieldti2, fieldti2(end-1:-1:1)], 2e3);
plot(0:0.2:2e3, iwi2)
plot(0:0.2:1e3, iwi2(1:5001))
607.8/0.2
fieldwh3a = dctI(fieldti2(1:3040))/sqrt(3039*pi);
fieldwh3a = dctI(fieldth3(1:2306))/sqrt(2305*pi);
fieldwi2a = dctI(fieldti2(1:3040))/sqrt(3039*pi);
fieldwi2b = dctI(fieldti2(3040:5001))/sqrt((5000 - 3040)*pi);
figure
plot(0:pi/607.8:pi/10, fieldwi2a(1:47))
length(0:pi/607.8:pi/10)
plot(0:pi/607.8:pi/10, fieldwi2a(1:61))
1e3-607.8
length(0:pi/392.2:pi/10)
hold on
plot(0:pi/392.2:pi/10, fieldwi2b(1:40))
0.78/0.048
evawi2a = dctI(evawi2(1:3040))/sqrt(3039*pi);
evawi2a = dctI(evati2(1:3040))/sqrt(3039*pi);
evawi2b = dctI(evati2(3040:5001))/sqrt((5000 - 3040)*pi);
figure
plot(0:pi/607.8:pi, fieldwi2a(1:601))
plot(0:pi/607.8:pi, fieldwi2a(1:600))
plot(0:pi/607.8:pi, fieldwi2a(1:602))
plot(0:pi/607.8:pi, fieldwi2a(1:610))
plot(0:pi/607.8:pi, fieldwi2a(1:611))
plot(0:pi/607.8:pi, fieldwi2a(1:609))
length(0:pi/607.8:pi)
plot(0:pi/607.8:pi, evawi2a(1:608))
hold on
plot(0:pi/392.2:pi, evawi2b(1:392))
plot(0:pi/392.2:pi, evawi2b(1:391))
plot(0:pi/392.2:pi, evawi2b(1:390))
length(0:pi/392.2:pi)
plot(0:pi/392.2:pi, evawi2b(1:393))
sum(psih3(:, end).*conj(psih3(:,end)))
sum(psib2(:, end).*conj(psib2(:,end)))
sum(psii2(:, end).*conj(psii2(:,end)))
J1i2
J1b2
figure
plot(0:0.2:1e3, psii2.*conj(psii2))
whos
[~, ~, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
psiEi2 = P240\psii2;
plot(0:0.2:1e3, psiEi2.*conj(psiEi2))
Vf
[~, E02402, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
E02402
xdomain240
Nx240
Vf
whos
Vf = @(x) 1 - 1/sqrt(1+x.^2);
[~, E02402, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
E0240
psiEi2 = P240\psii2;
plot(0:0.2:1e3, psiEi2.*conj(psiEi2))
max(P240(:, 1) - fi0240)
[~, E02402, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
figure
plot(x, V)
plot(x, Vf(x))
Vf = @(x) 1 - 1./sqrt(1+x.^2);
[~, E02402, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
E02402
psiEi2 = P240\psii2;
plot(0:0.2:1e3, psiEi2.*conj(psiEi2))
figure
npsii2 = sqnorm(psii2);
plot(0:0.2:1e3, npsii2)
E240(1:10)
E240(1:15)
E240(1:20)
E240(1:27)
E240(1:29)
figure
plot(x240, conj(P240(:,9).*P240(:,9)))
plot(x240, conj(P240(:,11).*P240(:,11)))
plot(x240, conj(P240(:,25).*P240(:,25)))
plot(x240, conj(P240(:,11).*P240(:,11)))
plot(x240, conj(P240(:,13).*P240(:,13)))
plot(x240, conj(P240(:,15).*P240(:,15)))
plot(x240, conj(P240(:,10).*P240(:,10)))
viewVPmiux(psii2, Vf, x240abs, fieldti2, x240, 0.01)
viewVPmiux(psii2, Vf, xabs240, fieldti2, x240, 0.01)
figure
plot(w(1:1001), evawi2)
plot(w(1:1001), evawi2(1:1001))
figure
plot(w(1:101), cww_con(1:101))
hold on
plot(w(1:101), fieldwi2(1:101))
figure
plot(0:0.2:1e3, iwi2(1:5001))