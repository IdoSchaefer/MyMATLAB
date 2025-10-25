t = 0:0.2:1e3;
plot(t, 0.1*exp(-(t-500).^2./(2*150^2)).*cos(0.06*(t-500)))
xlabel('$t$ (a.u.)', 'interpreter', 'latex')
ylabel('$\epsilon(t)$ (a.u.)', 'interpreter', 'latex')