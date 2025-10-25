dctfactor = 2000/sqrt(2000*pi)
[U1, mniter1, matvecs1] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.02*cos(0.06*t), [], fi0240, x240, 0:2000, 1e4, 7, 7, 1e-5);
[U2, mniter2, matvecs2] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.02*sin(0.06*t), [], fi0240, x240, 0:2000, 1e4, 7, 7, 1e-5);
ma1 = eva(U1, -x240./(1 + x240.^2).^(3/2), zeros(1, 2001));
ma1w = dctI(ma1)*dctfactor;
figure
plot(t, ma1)
plot(0:pi/2000:pi, log(ma1w(1:2001)))
ma2 = eva(U2, -x240./(1 + x240.^2).^(3/2), zeros(1, 2001));
ma2w = dctI(ma2)*dctfactor;
[U3, mniter3, matvecs3] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*0.02*sin(0.06*t).*(tanh((t-200)/70)+1), [], fi0240, x240, 0:2000, 1e4, 7, 7, 1e-5);
ma3 = eva(U3, -x240./(1 + x240.^2).^(3/2), zeros(1, 2001));
ma3w = dctI(ma3)*dctfactor;
w=0:pi/2000:pi;
w1=0:pi/2000:pi/0.2;
dctfactor1 = 2000/sqrt(1e4*pi)
t1=0:0.2:2000;
tf_t = (dctI(exp(-(w1-0.06).^2/0.01^2))/dctfactor1);
tfcos_t = dctI(exp(-(w1-0.06).^2/0.01^2).*cos(w1*2000))/dctfactor1;
a=tf_t(1)
b = tf_t(end)
field01 = fieldcos - (((fieldcos(1)*a - fieldcos(end)*b)*tf_t + (fieldcos(end)*a - fieldcos(1)*b)*tfcos_t)/(a^2-b^2));
figure
plot(t, field01)
plot(t1, field01)
[U4, mniter4, matvecs4] = SemiGlobalArnoldi_xp(K240, Vabs240, @(u,x,t) -xabs240*cosineIpln(field01, [0, 2000], t), [], fi0240, x240, 0:2000, 1e4, 7, 7, 1e-5);
ma4 = eva(U4, -x240./(1 + x240.^2).^(3/2), zeros(1, 2001));
ma4w = dctI(ma4)*dctfactor;