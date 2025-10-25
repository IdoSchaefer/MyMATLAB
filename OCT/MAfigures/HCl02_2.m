load HCl1
% J10 = 1.9223
% J1 = 1.6922e+001
[fieldtg, fieldwg, psig, evmiutg, evmiuwg, mnitercg, Jg, J1g, J2g, Jorthg, Jforbg] = guessresultsfalr(fi0, Vf, 1785, [-0.08 0.25],...
    xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2),...
    @(w) 100*rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 10000, 5, 5, 9, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfaloudrmiux(fi0, Vf, 1785, [-0.08 0.25],...
    xdomain, miuf3, P(:, 1:20), P(:, 21:32), @(w) rectanglefun(w, 0, 1.5e-2), @(w) 2500*rectanglefun(w, 0, 1.5e-2),...
    @(w) 100*rectanglefun(w, 2.5e-2, 2.7e-2), 0, @(n) n.^2, 1, 10000, 5, 5, 9, 1e-3);
w=0:pi/1e4:pi/5;
% The fieldw and evmiuw figure:
[AX, H1, H2]=plotyy(w, fieldw, w, evmiuw);
xlabel('$\omega$ (a.u.)', 'interpreter', 'latex')
% ylabel('$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
% ylabel('$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
% ylabel('mau')
% ylabel('$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
set(get(AX(2),'Ylabel'),'String','$\overline{\left<\hat\mu\right>}(\omega)$ (a.u.)', 'interpreter', 'latex')
set(get(AX(1),'Ylabel'),'String','$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
set(AX(2), 'FontSize', 16)
set(AX(2), 'Ycolor', 'b')
set(AX(2), 'Xlim', [0 0.027])
set(AX(2), 'Ylim', [-255 255])
set(get(AX(2),'YLabel'), 'Position', [0.028, 0, 0])
% The V(x) and miu(x) figure:
[AX, H1, H2] = plotyy(x, Vf(x), x, miuf3(x))
xlabel('$x$ (a.u.)', 'interpreter', 'latex')
set(get(AX(2),'Ylabel'),'String','$\mu(x)$ (a.u.)', 'interpreter', 'latex')
set(get(AX(1),'Ylabel'),'String','$V(x)$ (a.u.)', 'interpreter', 'latex')
set(AX(2), 'FontSize', 16)
set(AX(2), 'Ycolor', 'r')
set(AX(2), 'Xlim', xdomain)
set(AX(2), 'Ylim', [min(miuf3(x)) max(miuf3(x))])
% The fieldt, evmiut figure:
[AX, H1, H2]=plotyy(t, fieldt, t, evmiut);
set(get(AX(2),'Ylabel'),'String','$\left<\hat\mu\right>(t)$ (a.u.)', 'interpreter', 'latex')
set(get(AX(2),'YLabel'), 'Position', [0.028, 0, 0])
set(get(AX(1),'Ylabel'),'String','$\epsilon(t)$ (a.u.)', 'interpreter', 'latex')
set(AX(2), 'FontSize', 16)
set(AX(2), 'Ycolor', 'g')
set(AX(1), 'Ycolor', 'r')
set(AX(1), 'Ylim', [-0.166, 0.095])
set(AX(2), 'Ylim', [min(real(evmiut))-0.001, max(real(evmiut))+0.001])
xlabel('$t$ (a.u.)', 'interpreter', 'latex')
set(AX(2), 'Ytick', [-0.05, 0, 0.05 0.1])