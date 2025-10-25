% Spectra figure:
w=0:pi/1000:pi/0.2;
[AX, H1, H2]=plotyy(w, fieldw, w, evmiuw);
set(AX(1), 'FontSize', 16)
set(AX(1), 'Ycolor', 'r')
set(AX(1), 'Xlim', [0 0.6])
set(AX(2), 'Xlim', [0 0.6])
set(AX(2), 'FontSize', 16)
set(AX(2), 'Ycolor', 'b')
set(AX(2), 'Ylim', [-20 35])
set(AX(1), 'Ylim', [-20 35]/35*1.8)
set(AX(1), 'Ytick', [-1, -0.5, 0, 0.5, 1, 1.5])
set(AX(2), 'Ytick', [-20, -10, 0, 10, 20, 30])
xlabel('$\omega$ (a.u.)', 'interpreter', 'latex')
set(get(AX(2),'Ylabel'),'String','$\overline{\left<\hat\mu\right>}(\omega)$ (a.u.)', 'interpreter', 'latex')
set(get(AX(1),'Ylabel'),'String','$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
set(get(AX(1),'XLabel'), 'Position', [0, 0.28, 0])