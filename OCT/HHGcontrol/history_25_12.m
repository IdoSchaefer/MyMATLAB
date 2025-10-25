npsi3a = sqnorm(psi3a);
npsi6 = sqnorm(psi6);
npsi1b = sqnorm(psi1b);
npsi1b(end)
npsi3a(end)
npsi6(end)
npsicw1b(end)
npsicwe1b(end)
figure
plot(w(1:1001)/0.06, evaw1b(1:1001))
hold on
plot(w(1:1001)/0.06, evaw3a(1:1001))
npsi7a = sqnorm(psi7a);
npsi7a(end)
plot(w(1:1001)/0.06, evaw6(1:1001))
plot(w(1:1001)/0.06, evaw7a(1:1001))
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
ylabel('$\overline{\left<\hat\textbf{C}\right>}(\omega)$ (a.u.)','interpreter', 'latex')
clf
plot(w(1:1001)/0.06, evaw1b(1:1001))
hold on
plot(w(1:1001)/0.06, evaw3a(1:1001))
plot(w(1:1001)/0.06, evaw6(1:1001))
plot(w(1:1001)/0.06, evaw7a(1:1001))
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
ylabel('$\overline{\left<\hat\textbf{C}\right>}(\omega)$ (a.u.)','interpreter', 'latex')
figure
plot(t, fieldt1b)
xlabel('$t$ (a.u.)', 'interpreter', 'latex')
ylabel('$\epsilon(t)$ (a.u.)','interpreter', 'latex')
max(abs(fieldt7a))
hold on
sqrt(sum(fieldt1b.^2)/sum(cw_con.^2))
max(abs(fieldtcwe1b))
hold on
plot(t, fieldtcwe1b)
plot(t, fieldt3a)
clf
hold on
plot(t, fieldt3a)
hold on
plot(t, fieldt6)
hold on
plot(t, fieldt7a)
sum(fieldt1b.^2)*0.2
sum(fieldtcwe1b.^2)*0.2
sum(fieldtcw1b.^2)*0.2
sum(fieldtcw1b.^2)/sum(fieldtcwe1b.^2)
330.8/0.2
evaw3a_1 = dctI(fieldt3a(1:1655))/sqrt(pi*1655);
evaw3a_1 = dctI(fieldt3a(1:1655))/sqrt(pi*1654);
fieldt3a(1655)
t(1655)
5001-1654-1
evaw3a_2 = dctI(fieldt3a(1655:5001))/sqrt(pi*3346);
size(evaw3a_2)
w3a1 = 0:pi/330.8:pi/0.2;
size(w3a1)
1e3-330.8
w3a2 = 0:pi/669.2:pi/0.2;
size(w3a2)
figure
plot(w3a1, eva3a_1)
plot(w3a1, evaw3a_1)
hold on
plot(w3a2, evaw3a_2)
fieldt3a1 = dctI(fieldt3a(1:1655))/sqrt(pi*1654);
clear fieldt3a1;
fieldw3a1 = dctI(fieldt3a(1:1655))/sqrt(pi*1654);
fieldw3a2 = dctI(fieldt3a(1655:5001))/sqrt(pi*3346);
clf
plot(w3a1, fieldw3a1)
hold on
plot(w3a2, fieldw3a2)
plot(w, fieldw3a/1e3)
clf
plot(w3a1/0.06, fieldw3a1)
hold on
plot(w3a2/0.06, fieldw3a2)
plot(w3a2/0.06, fieldw3a2*max(abs(fieldw3a2))/max(abs(fieldw3a1)))
plot(w3a2/0.06, fieldw3a2*max(abs(fieldw3a1))/max(abs(fieldw3a2)))
clear evaw3a_1 evaw3a_2
evaw3a_1 = dctI(evat3a(1:1655))/sqrt(pi*1654);
clear evaw3a_1 evaw3a_2
evaw3a1 = dctI(evat3a(1:1655))/sqrt(pi*1654);
evaw3a2 = dctI(evat3a(1655:5001))/sqrt(pi*3346);
figure
plot(w3a1/0.06, evaw3a1)
hold on
plot(w3a2/0.06, evaw3a2)
clf
plot(w3a1/0.06, evaw3a1)
hold on
plot(w3a1/0.06, evaw3a2)
plot(w3a2/0.06, evaw3a2)
fieldw1b1 = dctI(fieldt3a(1:1655))/sqrt(pi*1654
348.2/0.2
fieldw1b1 = dctI(fieldt1b(1:1742))/sqrt(pi*1741);
fieldw1b2 = dctI(fieldt1b(1742:5001))/sqrt(pi*
5001-1741-1
fieldw1b2 = dctI(fieldt1b(1742:5001))/sqrt(pi*3259);
figure
w1b1 = 0:pi/348.2:pi/0.2;
w1b2 = 0:pi/651.8:pi/0.2;
size(w1b1)
size(w1b2)
1e3-348.2
figure
plot(w1b1/0.06, fieldw1b1)
hold on
plot(w1b2/0.06, fieldw1b2)
figure
evaw1b1 = dctI(evat1b(1:1742))/sqrt(pi*1741);
evaw1b2 = dctI(evat1b(1742:5001))/sqrt(pi*3259);
figure
plot(w1b1/0.06, evaw1b1)
hold on
plot(w1b2/0.06, evaw1b2)
save fields_13_15_17
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
ylabel('$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
ylabel('$\overline{\left<\hat\textbf{C}\right>}(\omega)$ (a.u.)','interpreter', 'latex')
hold on
plot(348.2, [-0.1 0.1])
plot([348.2, 348.2], [-0.1 0.1])
ylabel('$\overline{\left<\hat\textbf{C}\right>}(\omega)/\Delta t$ (a.u.)','interpreter', 'latex')
ylabel('$\bar\epsilon(\omega)/\Delta t$ (a.u.)','interpreter', 'latex')
plot(w1b2/0.06, fieldw1b2*348.2)
figure
plot(w1b2/0.06, fieldw1b2*348.2)
hold on
clf
1e3-348.2
plot(w1b1/0.06, fieldw1b1*348.2)
hold on
plot(w1b2/0.06, fieldw1b2*651.8)
xlabel('$\omega/\omega_0$', 'interpreter', 'latex')
ylabel('$\bar\epsilon(\omega)$ (a.u.)','interpreter', 'latex')