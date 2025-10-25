figure
plot(w(1:1001)/0.06, evaw1b(1:1001))
hold on
plot(w(1:1001)/0.06, evawcw1b(1:1001))
plot(w(1:1001)/0.06, evawcwe1b(1:1001))
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
ylabel('$\overline{\left<\hat\textbf{C}\right>}(\omega)$ (a.u.)','interpreter', 'latex')