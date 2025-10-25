figure
yyaxis left
plot(w(1:1001)/0.06, fieldw1b(1:1001))
yyaxis right
plot(w(1:1001)/0.06, evaw1b(1:1001))
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
yyaxis left
ylabel('$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
yyaxis right
ylabel('$\overline{\left<\hat\textbf{C}\right>}(\omega)$ (a.u.)','interpreter', 'latex')
