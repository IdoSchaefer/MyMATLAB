load TLS2pi
whos
nu_cut = 26.64/2.5
stiffness = 0.1155*2.5
plot(0:0.5:50, filterEcut(0:0.5:50))
ftarget
Hoperations_exp
muxz_exp
figure
plot(0:0.5:50, filterEcut(0:0.5:50).^2)
22.81/2.5
figure
plot(0:0.01:1, fieldtf1)
plot(0:0.01:1, fieldtf)
plot(0:0.5:50, dctIcos2pi)
H0
[Hoperations_exp1, fcouplingOp_exp1] = Hmats2Hops(Sz, muxz_exp)
Sz = [1/2 0; 0 -1/2]
muxz_exp
[Hoperations_exp1, fcouplingOp_exp1] = Hmats2Hops(Sz, muxz_exp)
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
ftarget
ftarget1 = get_fchi_state([1;0])
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
figure
plot(0:0.01:1, fieldtf1)
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
figure
plot(0:0.01:1, fieldtf1)
[fieldtfSz, fieldnufSz, psifSz, relEfSz, convfSz, niterfSz, mallnitercfSz, JtermsfSz, maxgradfSz, alphafSz, invHessfSz] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
Jtermsf1
1-Jtermsf1.Jmax
JtermsfSz
1-JtermsfSz.Jmax
figure
plot(0:0.01:1, fieldtfSz)
optionspi
[fieldtfSz, fieldnufSz, psifSz, relEfSz, convfSz, niterfSz, mallnitercfSz, JtermsfSz, maxgradfSz, alphafSz, invHessfSz] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
figure
plot(0:0.01:1, fieldtfSz)
hold on
plot(0:0.01:1, fieldtf1)
plot(0:0.01:1, -fliplr(fieldtf1))
JtermsfSz
Jtermsf1
figure
plot(0:0.01:1, -fieldtfSz)
hold on
plot(0:0.01:1, fieldtf1)
figure
plot(0:0.01:1, fieldtf1+fliplr(fieldtfSz))
convfz(end)
convfSz(end)
convfSz(end)-convf1
convfSz(end)-convf1(end)
figure
plot(0:0.01:1, -fieldtfSz)
hold on
plot(0:0.01:1, -fieldtfSz*2)
clf
plot(0:0.01:1, -fieldtfSz*2)
hold on
plot(0:0.01:1, fieldtf1*2)
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
legend('$\hat{\textbf{H}}_r(\tau)=\frac{1}{2}\sigma_z + \frac{1}{2}\epsilon(\tau)\left(\sigma_x + \frac{\sigma_z}{\sqrt{2}}\right)$','$\hat{\textbf{H}}_r(\tau)=-\frac{1}{2}\sigma_z - \frac{1}{2}\epsilon(\tau)\left(\sigma_x + \frac{\sigma_z}{\sqrt{2}}\right)$','interpreter','latex')
figure
plot(0:0.01:1, 2*(fieldtf1+fliplr(fieldtfSz)))
xlabel('$\frac{\tau}{2\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon_2(\tau)-\epsilon_1(2\pi-\tau)$', 'interpreter', 'latex')
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], 2/3*dctIcos2pi(1:41), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
optionspi
dctIcos3pi = dctI(cos(0:pi/50:3*pi))*sqrt(pi/100);
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos3pi(1:41), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
figure
plot(0:0.5:50, dctIcos2pi)
plot(0:1/3:50, dctIcos3pi)
dctIcos3pi(abs(dctIcos3pi)<1e-8)= 0;
plot(0:1/3:50, dctIcos3pi)
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos3pi(1:41), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos3pi(1:61), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
figure
[fieldtf1, fieldnuf1, psif1, relEf1, convf1, niterf1, mallnitercf1, Jtermsf1, maxgradf1, alphaf1, invHessf1] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos2pi(1:41), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
[fieldtf3pi, fieldnuf3pi, psif3pi, relEf3pi, convf3pi, niterf3pi, mallnitercf3pi, Jtermsf3pi, maxgradf3pi, alphaf3pi, invHessf3pi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos3pi(1:61), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
Jtermsf3pi
1-Jtermsf3pi.Jmax
figure
plot(0:1/3:50, fieldtf3pi)
plot(0:1/150:1, fieldtf3pi)
xlabel('$\frac{\tau}{3\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
plot(0:1/150:1, 3*fieldtf3pi)
xlabel('$\frac{\tau}{3\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
dctIcos3pi = dctI(cos(0:pi/50:3*pi))*sqrt(pi/150);
[fieldtf3pi, fieldnuf3pi, psif3pi, relEf3pi, convf3pi, niterf3pi, mallnitercf3pi, Jtermsf3pi, maxgradf3pi, alphaf3pi, invHessf3pi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos3pi(1:61), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
dctIcos3pi(abs(dctIcos3pi)<1e-8)= 0;
[fieldtf3pi, fieldnuf3pi, psif3pi, relEf3pi, convf3pi, niterf3pi, mallnitercf3pi, Jtermsf3pi, maxgradf3pi, alphaf3pi, invHessf3pi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos3pi(1:61), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
Jtermsf3pi
1-Jtermsf3pi.Jmax
plot(0:1/150:1, 3*fieldtf3pi)
xlabel('$\frac{\tau}{3\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
dctIcos4pi = dctI(cos(0:pi/50:4*pi))*sqrt(pi/200);
dctIcos4pi(abs(dctIcos4pi)<1e-8)= 0;
[fieldtf4pi, fieldnuf4pi, psif4pi, relEf4pi, convf4pi, niterf4pi, mallnitercf4pi, Jtermsf4pi, maxgradf4pi, alphaf4pi, invHessf4pi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos4pi(1:61), @(nu) 200*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 4*pi, pi/50, 7, 7, 1e-7);
[fieldtf3pi, fieldnuf3pi, psif3pi, relEf3pi, convf3pi, niterf3pi, mallnitercf3pi, Jtermsf3pi, maxgradf3pi, alphaf3pi, invHessf3pi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -dctIcos3pi(1:61), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
plot(0:1/150:1, 3*fieldtf3pi)
plot(0:1/150:1, -3*fieldtf3pi)
xlabel('$\frac{\tau}{3\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
[fieldtf3pi21, fieldnuf3pi21, psif3pi21, relEf3pi21, convf3pi21, niterf3pi21, mallnitercf3pi21, Jtermsf3pi21, maxgradf3pi21, alphaf3pi21, invHessf3pi21] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos3pi(1:61), @(nu) 150*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3*pi, pi/50, 7, 7, 1e-7);
figure
plot(0:1/150:1, 3*fieldtf3pi21)
hold on
plot(0:1/150:1, -3*fieldtf3pi)
plot(0:1/150:1, 3*fliplr(fieldtf3pi21))
figure
plot(0:1/150:1, 3*(fliplr(fieldtf3pi21)+ fieldtf3pi21)
plot(0:1/150:1, 3*(fliplr(fieldtf3pi21)+ fieldtf3pi21))
plot(0:1/150:1, 3*(fliplr(fieldtf3pi21)- fieldtf3pi21))
figure
plot(0:1/150:1, 3*(fieldtf3pi + fieldtf3pi21))
xlabel('$\frac{\tau}{3\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
[fieldtf4pi, fieldnuf4pi, psif4pi, relEf4pi, convf4pi, niterf4pi, mallnitercf4pi, Jtermsf4pi, maxgradf4pi, alphaf4pi, invHessf4pi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], dctIcos4pi(1:81), @(nu) 200*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 4*pi, pi/50, 7, 7, 1e-7);
figure
plot(0:1/200:1, -4*fieldtf4pi)
xlabel('$\frac{\tau}{4\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
[fieldtf4pi21, fieldnuf4pi21, psif4pi21, relEf4pi21, convf4pi21, niterf4pi21, mallnitercf4pi21, Jtermsf4pi21, maxgradf4pi21, alphaf4pi21, invHessf4pi21] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos3pi(1:81), @(nu) 200*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 4*pi, pi/50, 7, 7, 1e-7);
[fieldtf4pi21, fieldnuf4pi21, psif4pi21, relEf4pi21, convf4pi21, niterf4pi21, mallnitercf4pi21, Jtermsf4pi21, maxgradf4pi21, alphaf4pi21, invHessf4pi21] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], dctIcos4pi(1:81), @(nu) 200*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 4*pi, pi/50, 7, 7, 1e-7);
figure
plot(0:1/200:1, 4*fieldtf4pi21)
[fieldtf4pi21, fieldnuf4pi21, psif4pi21, relEf4pi21, convf4pi21, niterf4pi21, mallnitercf4pi21, Jtermsf4pi21, maxgradf4pi21, alphaf4pi21, invHessf4pi21] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], -dctIcos4pi(1:81), @(nu) 200*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 4*pi, pi/50, 7, 7, 1e-7);
figure
plot(0:1/200:1, 4*fieldtf4pi21)
figure
plot(0:1/200:1, 4*(fieldtf4pi21-fieldtf4pi))
plot(0:1/200:1, 4*(fieldtf4pi21+fieldtf4pi))
1-Jtermsf4pi.Jmax
optionspig = optionspi
optionspig.maxNiter = 0
[fieldtf1gi, fieldnuf1gi, psif1gi, relEf1gi, convf1gi, niterf1gi, mallnitercf1gi, Jtermsf1gi, maxgradf1gi, alphaf1gi, invHessf1gi] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], fieldnuf1(end:-1:1), @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
[fieldtf1gi, fieldnuf1gi, psif1gi, relEf1gi, convf1gi, niterf1gi, mallnitercf1gi, Jtermsf1gi, maxgradf1gi, alphaf1gi, invHessf1gi] = OClimf_qn([1;0], ftarget, Hoperations_exp, 1, fcouplingOp_exp, [-6 6], -fieldnufSz, @(nu) 100*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspig, 2*pi, 2*pi/1e2, 7, 7, 1e-7);
Jtermsf1gi
Jtermsf1gi.Jmax - Jtermsf1.Jmax
convf1(end)-convf1gi
save TLS2pi
xlabel('$\frac{\tau}{4\pi}$', 'interpreter', 'latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
[fieldtfpi, fieldnufpi, psifpi, relEfpi, convfpi, niterfpi, mallnitercfpi, Jtermsfpi, maxgradfpi, alphafpi, invHessfpi] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*sin(nu), @(nu) 50*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, pi, pi/50, 7, 7, 1e-7);
figure
plot(0:pi/50:1.2*pi, fieldtfpi)
size(fieldtfpi)
[fieldtfpi2dt, fieldnufpi2dt, psifpi2dt, relEfpi2dt, convfpi2dt, niterfpi2dt, mallnitercfpi2dt, Jtermsfpi2dt, maxgradfpi2dt, alphafpi2dt, invHessfpi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*sin(nu), @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt.Jmax
[fieldtfpi2dt, fieldnufpi2dt, psifpi2dt, relEfpi2dt, convfpi2dt, niterfpi2dt, mallnitercfpi2dt, Jtermsfpi2dt, maxgradfpi2dt, alphafpi2dt, invHessfpi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.05*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*sin(nu), @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt.Jmax
[fieldtfpi2dt, fieldnufpi2dt, psifpi2dt, relEfpi2dt, convfpi2dt, niterfpi2dt, mallnitercfpi2dt, Jtermsfpi2dt, maxgradfpi2dt, alphafpi2dt, invHessfpi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*sin(nu), @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
figure
plot(0:pi/50:pi, fieldtfpi)
figure
plot(0:pi/50:1.2*pi, fieldtfpi2dt)
[fieldtfpi2dt1, fieldnufpi2dt1, psifpi2dt1, relEfpi2dt1, convfpi2dt1, niterfpi2dt1, mallnitercfpi2dt1, Jtermsfpi2dt1, maxgradfpi2dt1, alphafpi2dt1, invHessfpi2dt1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*sin(nu), @(nu) 6*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt1.Jmax
Jtermsfpi2dt.Jmax
figure
plot(0:pi/50:1.2*pi, fieldtfpi2dt1)
plot(0:pi/50:1.2*pi, -fieldtfpi2dt1)
plot((0:1/50:1.2)/2, -fieldtfpi2dt1)
plot((0:1/50:1.2)/2, -fieldtfpi2dt)
[fieldtfpi2dt2, fieldnufpi2dt2, psifpi2dt2, relEfpi2dt2, convfpi2dt2, niterfpi2dt2, mallnitercfpi2dt2, Jtermsfpi2dt2, maxgradfpi2dt2, alphafpi2dt2, invHessfpi2dt2] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*sin(nu), @(nu) 600*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt1.Jmax
Jtermsfpi2dt2.Jmax
[fieldtfpi2dt2, fieldnufpi2dt2, psifpi2dt2, relEfpi2dt2, convfpi2dt2, niterfpi2dt2, mallnitercfpi2dt2, Jtermsfpi2dt2, maxgradfpi2dt2, alphafpi2dt2, invHessfpi2dt2] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 600*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt2.Jmax
figure
plot((0:1/50:1.2)/2, -fieldtfpi2dt2)
1-Jtermsfpi2dt2.Jmax
[fieldtfpi2dt2, fieldnufpi2dt2, psifpi2dt2, relEfpi2dt2, convfpi2dt2, niterfpi2dt2, mallnitercfpi2dt2, Jtermsfpi2dt2, maxgradfpi2dt2, alphafpi2dt2, invHessfpi2dt2] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.01)), @(nu) 600*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
[fieldtfpi2dt2, fieldnufpi2dt2, psifpi2dt2, relEfpi2dt2, convfpi2dt2, niterfpi2dt2, mallnitercfpi2dt2, Jtermsfpi2dt2, maxgradfpi2dt2, alphafpi2dt2, invHessfpi2dt2] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25^2)), @(nu) 600*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt2.Jmax
[fieldtfpi2dt2, fieldnufpi2dt2, psifpi2dt2, relEfpi2dt2, convfpi2dt2, niterfpi2dt2, mallnitercfpi2dt2, Jtermsfpi2dt2, maxgradfpi2dt2, alphafpi2dt2, invHessfpi2dt2] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.5^2)), @(nu) 600*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt2.Jmax
figure
Jtermsfpi2dt2.Jenergy
plot((0:1/50:1.2)/2, -fieldtfpi2dt2)
[fieldtfpi2dt1, fieldnufpi2dt1, psifpi2dt1, relEfpi2dt1, convfpi2dt1, niterfpi2dt1, mallnitercfpi2dt1, Jtermsfpi2dt1, maxgradfpi2dt1, alphafpi2dt1, invHessfpi2dt1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 6*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt1.Jenergy
Jtermsfpi2dt1.Jenergymax
Jtermsfpi2dt1.Jmax
[fieldtfpi2dt, fieldnufpi2dt, psifpi2dt, relEfpi2dt, convfpi2dt, niterfpi2dt, mallnitercfpi2dt, Jtermsfpi2dt, maxgradfpi2dt, alphafpi2dt, invHessfpi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt.Jmax
[fieldtf3pi2dt, fieldnuf3pi2dt, psif3pi2dt, relEf3pi2dt, convf3pi2dt, niterf3pi2dt, mallnitercf3pi2dt, Jtermsf3pi2dt, maxgradf3pi2dt, alphaf3pi2dt, invHessf3pi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 160*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 3.2*pi, pi/50, 7, 7, 1e-7);
Jtermsf3pi2dt.Jmax
1-Jtermsf3pi2dt.Jmax
figure
plot((0:1/50:3.2)/2, -fieldtf3pi2dt)
figure
plot((0:1/50:3)/2, -fieldtf3pi)
whos
[fieldtf10pi2dt, fieldnuf10pi2dt, psif10pi2dt, relEf10pi2dt, convf10pi2dt, niterf10pi2dt, mallnitercf10pi2dt, Jtermsf10pi2dt, maxgradf10pi2dt, alphaf10pi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 510*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 10.2*pi, pi/50, 7, 7, 1e-7);
[fieldtf10pi2dt, fieldnuf10pi2dt, psif10pi2dt, relEf10pi2dt, convf10pi2dt, niterf10pi2dt, mallnitercf10pi2dt, Jtermsf10pi2dt, maxgradf10pi2dt, alphaf10pi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.01*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 510*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 10.2*pi, pi/50, 7, 7, 1e-7);
[fieldtf10pi2dt, fieldnuf10pi2dt, psif10pi2dt, relEf10pi2dt, convf10pi2dt, niterf10pi2dt, mallnitercf10pi2dt, Jtermsf10pi2dt, maxgradf10pi2dt, alphaf10pi2dt] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) (1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 510*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 10.2*pi, pi/50, 7, 7, 1e-7);
Jtermsf10pi2dt.Jmax
1-Jtermsf10pi2dt.Jmax
figure
plot((0:1/50:10.2)/2, -fieldtf10pi2dt)
[fieldtfpi2dtg, fieldnufpi2dtg, ~, ~, ~, ~, ~, Jtermsfpi2dtg] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/(2*0.25)), @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspig, 1.2*pi, pi/50, 7, 7, 1e-7);
figure
plot((0:1/50:1.2)/2, -fieldtfpi2dtg)
figure
plot(0:22, fieldnufpi2dtg(1:23))
optionspig
figure
plot(0:22, 0.1*(1-tanh(stiffness*((0:22)-nu_cut))).*(1-heaviside((0:22)-20.1)).*exp(-((0:22)-1).^2/(2*0.25)))
[fieldtfpi2dtg1, fieldnufpi2dtg1, ~, ~, ~, ~, ~, Jtermsfpi2dtg1] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/2).*nu, @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspig, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dtg1
Jtermsfpi2dtg
figure
plot((0:1/50:1.2)/2, -fieldtfpi2dtg1)
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], @(nu) 0.1*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)).*exp(-(nu-1).^2/2).*nu, @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt3
dctIsin_pi = dctI(sin((0:pi/50:1.2*pi)-0.6*pi))*sqrt(pi/50);
figure
plot((0:pi/50:1.2*pi), sin((0:pi/50:1.2*pi)-0.6*pi))
clear dctIsin_pi
dctIsin_pi = dctI(cos((0:pi/50:1.2*pi)-0.6*pi))*sqrt(pi/50);
plot((0:pi/50:1.2*pi), cos((0:pi/50:1.2*pi)-0.6*pi))
figure
plot(0:22, dctIsin_pi(1:23))
fieldnupi2dtg_sin = dctIsin_pi.*(1-tanh(stiffness*((0:50)-nu_cut))).*(1-heaviside((0:50)-20.1));
size(dctIsin_pi)
fieldnupi2dtg_sin = dctIsin_pi.*(1-tanh(stiffness*((0:5/6:50)-nu_cut))).*(1-heaviside((0:5/6:50)-20.1));
plot(0:5/6:50, dctIsin_pi)
figure
plot(0:5/6:50, fieldnupi2dtg_sin)
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], fieldnupi2dtg_sin, @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt3
Jtermsfpi2dt
figure
plot(0:pi/50:1.2*pi, -fieldtfpi2dt3)
plot((0:1/50:1.2)/2, -fieldtfpi2dt3)
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -fieldnupi2dtg_sin, @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt
Jtermsfpi2dt3
figure
plot((0:1/50:1.2)/2, -fieldtfpi2dt3)
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -0.1*fieldnupi2dtg_sin, @(nu) 60*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -fieldnupi2dtg_sin, @(nu) 6*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -0.1*fieldnupi2dtg_sin, @(nu) 6*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt3
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -0.1*fieldnupi2dtg_sin, @(nu) 20*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt3
[fieldtfpi2dt3, fieldnufpi2dt3, psifpi2dt3, relEfpi2dt3, convfpi2dt3, niterfpi2dt3, mallnitercfpi2dt3, Jtermsfpi2dt3, maxgradfpi2dt3, alphafpi2dt3, invHessfpi2dt3] = OClimf_qn([0;1], ftarget1, Hoperations_exp1, 1, fcouplingOp_exp1, [-6 6], -fieldnupi2dtg_sin, @(nu) 20*0.5*(1-tanh(stiffness*(nu-nu_cut))).*(1-heaviside(nu-20.1)), optionspi, 1.2*pi, pi/50, 7, 7, 1e-7);
Jtermsfpi2dt3
figure
plot((0:1/50:1.2)/2, -fieldtfpi2dt3)
Jtermsfpi2dt
Jtermsfpi2dt1
Jtermsfpi2dt2
1-Jtermsfpi2dt.Jmax
1-Jtermsfpi2dt1.Jmax
1-Jtermsfpi2dt2.Jmax
figure
plot((0:1/50:1.2)/2, -fieldtfpi2dt)
clf
plot((0:1/50:1.2)/2, -fieldtfpi2dt2)
hold on
plot((0:1/50:1.2)/2, -fieldtfpi2dt1)
plot((0:1/50:1.2)/2, -fieldtfpi2dt)
plot((0:1/50:1.2)/2, -fieldtfpi2dt1)
xlabel('$\frac{\omega_0 t}{2\pi}$', 'interpreter', 'latex')
ylabel('drive amplitude', 'interpreter', 'latex')
legend
save TLS2pi
figure
plot((0:1/50:3.2)/2, -fieldtf3pi2dt)
xlabel('$\frac{\omega_0 t}{2\pi}$', 'interpreter', 'latex')
ylabel('drive amplitude', 'interpreter', 'latex')
figure
plot((0:1/50:10.2)/2, -fieldtf10pi2dt)
xlabel('$\frac{\omega_0 t}{2\pi}$', 'interpreter', 'latex')
ylabel('drive amplitude', 'interpreter', 'latex')
save drives_exp_limf_1_3_10 fieldtfpi2dt fieldtfpi2dt1 fieldtfpi2dt2 fieldtf3pi2dt fieldtf10pi2dt fieldnufpi2dt fieldnufpi2dt1 fieldnufpi2dt2 fieldnuf3pi2dt fieldnuf10pi2dt Jtermsfpi2dt Jtermsfpi2dt1 Jtermsfpi2dt2 Jtermsf3pi2dt Jtermsf10pi2dt
1-Jtermsfpi2dt.Jmax
1-Jtermsfpi2dt1.Jmax
1-Jtermsfpi2dt2.Jmax
1-Jtermsf3pi2dt.Jmax
1-Jtermsf10pi2dt.Jmax
nucut
nu_cut
stiffness
nu_grid = 0:0.1:20;
figure
plot(nu_grid, 0.5*(1-tanh(stiffness*(nu_grid-nu_cut)))
plot(nu_grid, 0.5*(1-tanh(stiffness*(nu_grid-nu_cut))))
xlabel('$\frac{\omega}{\omega_0}$', 'interpreter', 'latex')
ylabel('\phi(\omega)', 'interpreter', 'latex')
ylabel('$\phi(\omega)$', 'interpreter', 'latex')
1-JtermsfSz.Jmax
1-Jtermsf3pi.Jmax
1-Jtermsf4pi.Jmax
optionspi
optionspi.f_termination