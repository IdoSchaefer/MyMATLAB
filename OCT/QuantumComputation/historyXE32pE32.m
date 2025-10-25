load('qubitH.mat')
whos
x32
[fi032, E032, x32, E32, P32, H32] = gsV(@(x) 26.45*pi*(1-cos(x)), [-pi pi], 32, m
xdlength=2*pi;
dx = xdlength/32
p = (-pi/dx:(2*pi/xdlength):(pi/dx - 2*pi/xdlength)).';
K = p.^2/(2*m);
m=1/(3.2*pi)
clear p
p32 = (-pi/dx:(2*pi/xdlength):(pi/dx - 2*pi/xdlength)).';
clear dx
dx32 = xdlength/32
p32 = (-pi/dx:(2*pi/xdlength):(pi/dx - 2*pi/xdlength)).';
p32 = (-pi/dx32:(2*pi/xdlength):(pi/dx32 - 2*pi/xdlength)).';
px32 = p2x(diag(p));
px32 = p2x(diag(p32));
pE32 = P32\px32*P32;
XE32(1:7, 1:7)
pE32(1:7, 1:7)
clear pE
XE32(1:7, 1:7)./abs(pE(1:7, 1:7)
XE32(1:7, 1:7)./abs(pE(1:7, 1:7))
XE32(1:7, 1:7)./abs(pE32(1:7, 1:7))
px(1:7,1:7)
px32(1:7,1:7)
(XE32*pE32-pE32*XE32)
px32*x32
px32/1i*x32
px32/1i*[x32(17:end);x32(1:16)]
edit p2x
p32
ifft([p32(17:end);p32(1:16)]*fft(x32))
ifft([p32(17:end);p32(1:16)].*fft(x32))
ifft([p32(17:end);p32(1:16)].^2.*fft(x32))
XE32(1:7, 1:7)./sqrt(m*4.6*2*pi)
XE32(1:7, 1:7)/sqrt(m*4.6*2*pi)
XE32(1:7, 1:7)/sqrt(m*4.6*2*pi)*sqrt(2)
XE32(1:7, 1:7)*sqrt(m*4.6*2*pi)*sqrt(2)
pE32(1:7, 1:7)/sqrt(m*4.6*2*pi)*sqrt(2)
sqrt(2:7))
sqrt(2:7)
figure
plot(x32, conj(P32(:,1).*P32(:,1))
plot(x32, conj(P32(:,1)).*P32(:,1))
plot(x32, conj(P32(:,2)).*P32(:,2))
XE32 = P32'*diag(x32)*P32;
XE32(1:7, 1:7)
XE32(1:7, 1:7)*sqrt(m*4.6*2*pi)*sqrt(2)
pE32 = P32'*px32*P32;
pE32(1:7, 1:7)/sqrt(m*4.6*2*pi)*sqrt(2)
plot(x32, conj(P32(:,3)).*P32(:,3))
plot(x32, conj(P32(:,4)).*P32(:,4))
plot(x32, conj(P32(:,5)).*P32(:,5))
P32ext = [P32(1, :)/sqrt(2); P32(2:32, :); P32(1, :)/sqrt(2)];
sqnorm(P32ext)
x32ext = [x32(1:32); -x32(1)]
XE32ext = P32ext'*diag(x32ext)*P32ext
XE32(1:7, 1:7)*sqrt(m*4.6*2*pi)*sqrt(2)
XE32ext(1:7, 1:7)*sqrt(m*4.6*2*pi)*sqrt(2)
whos
save qubitH
xlabel('$t$ (ns)', 'interpreter', 'latex')
%-- 09/07/2020 9:21 --%
load('qubitH.mat')
load('unitary3LS.mat')
figure
plot(0:0.01:10, fieldu1)
xlabel('$t$ (ns)', 'interpreter', 'latex')
ylabel('$\epsilon(t)$',  'interpreter', 'latex')
XE32ext(1:7, 1:7)*sqrt(m*4.6*2*pi)*sqrt(2)
sqrt(2:7)
pE32(1:7, 1:7)/sqrt(m*4.6*2*pi)*sqrt(2)