% Spectra figure:
w=0:pi/2000:pi/0.2;
[AX, H1, H2]=plotyy(w, fieldw, w, evmiuw);
set(AX(1), 'FontSize', 16)
set(AX(1), 'Ycolor', 'r')
set(AX(1), 'Xlim', [0 0.7])
set(AX(2), 'Xlim', [0 0.7])
set(AX(2), 'FontSize', 16)
set(AX(2), 'Ycolor', 'b')
set(AX(2), 'Ylim', [-52 57])
set(AX(1), 'Ylim', [-52 57]/57)
set(AX(2), 'Ytick', [-50, -25, 0, 25, 50])
xlabel('$\omega$ (a.u.)', 'interpreter', 'latex')
set(get(AX(2),'Ylabel'),'String','$\overline{\left<\hat\mu\right>}(\omega)$ (a.u.)', 'interpreter', 'latex')
set(get(AX(1),'Ylabel'),'String','$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
set(get(AX(1),'XLabel'), 'Position', [0, 0.28, 0])
% s(x), gamma(x) figure:
figure
[AX, H1, H2]=plotyy(x, 0.5*(tanh(x+35)-tanh(x-35)), x, 1e-3*(-0.5*(tanh(x+35)-tanh(x-35))+1));
xlabel('$x$ (a.u.)', 'interpreter', 'latex')
%set(AX(1), 'Ycolor', 'r')
set(AX(1), 'Ycolor', 'b')
set(AX(2), 'Ycolor', 'g')
set(get(AX(1),'Ylabel'),'String','$s(x)$','interpreter', 'latex')
set(get(AX(2),'Ylabel'),'String','$\gamma(x)$', 'interpreter', 'latex')
set(AX(1), 'FontSize', 16)
set(AX(2), 'FontSize', 16)