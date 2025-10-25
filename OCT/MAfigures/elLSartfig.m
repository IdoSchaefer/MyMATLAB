w=0:pi/100:20*pi;
% The fieldw and evmiuw figure:
[AX, H1, H2]=plotyy(w, fieldw, w, evmiuw);
xlabel('$\omega$ (a.u.)', 'interpreter', 'latex')
set(get(AX(2),'Ylabel'),'String','$\overline{\left<\hat\mu\right>}(\omega)$ (a.u.)', 'interpreter', 'latex')
set(get(AX(1),'Ylabel'),'String','$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
set(AX(1), 'FontSize', 16)
set(AX(1), 'Ycolor', 'r')
set(AX(1), 'Xlim', [0 12])
set(AX(1), 'Ylim', [-17/12*4.2 4.2])
set(AX(2), 'FontSize', 16)
set(AX(2), 'Ycolor', 'b')
set(AX(2), 'Xlim', [0 12])
set(AX(2), 'Ylim', [-17 12])
set(get(AX(2),'YLabel'), 'Position', [0.028, 0, 0])
set(AX(1), 'Ytick', [-4, -2, 0, 2, 4])
set(AX(2), 'Ytick', [-15, -10, -5, 0, 5, 10])