[fi0, E0, phi, E, P, H] = gsV(@(phi)2.5*(1-cos(phi)), [-pi pi], Nphi, 2.5);
[fi0, E0, phi, E, P, H] = gsV(@(phi)2.5*(1-cos(phi)), [-pi pi], 10, 2.5);
E
[fi016, E016, phi16, E16, P16, H16] = gsV(@(phi)2.5*(1-cos(phi)), [-pi pi], 16, 2.5);
[E E16(1:10)]
[fi032, E032, phi32, E32, P32, H32] = gsV(@(phi)2.5*(1-cos(phi)), [-pi pi], 32, 2.5);
[E E16(1:10) E32(1:10)]
figure
plot(phi32, P32(:,1))
plot(phi32, P32(:,1).*conj(P32(:,1)))
plot(phi32, P32(:,2).*conj(P32(:,2)))
plot(phi32, P32(:,3).*conj(P32(:,3)))
plot(phi32, P32(:,4).*conj(P32(:,4)))
plot(phi32, P32(:,5).*conj(P32(:,5)))
plot(phi32, P32(:,6).*conj(P32(:,6)))
plot(phi32, P32(:,7).*conj(P32(:,7)))
plot(phi32, P32(:,8).*conj(P32(:,8)))
plot(phi32, P32(:,6).*conj(P32(:,6)))
hold on
plot(phi32, P16(:,6).*conj(P16(:,6)))
plot(phi16, P16(:,6).*conj(P16(:,6)))
plot(phi16, P16(:,6).*conj(P16(:,6))/sqrt(2))
plot(phi16, P16(:,6).*conj(P16(:,6))/2)
figure
plot(phi16, P16(:,6).*conj(P16(:,6))/2-P32(1:2:end,6).*conj(P32(1:2:end,6)))
(E16(1:10)-E32(1:10))./E32(1:10)
(E-E32(1:10))./E32(1:10)
n=1:7;
0.5*n.*(n-1)
n=(1:7).';
n-0.5*n.*(n-1)/20
n-0.5*n.*(n-1)/20-0.5
[n-0.5*n.*(n-1)/20-0.5 E32(1:7)]
EDuffing = n-0.5*n.*(n-1)/20-0.5
[EDuffing(2:7) - Eduffing(1:6), E32(2:7)-E32(1:6)]
[EDuffing(2:7) - EDuffing(1:6), E32(2:7)-E32(1:6)]
n=(0:6).';
EDuffing = n-0.5*n.*(n-1)/20+0.5
[EDuffing(2:7) - EDuffing(1:6), E32(2:7)-E32(1:6)]
[(EDuffing(2:7) - EDuffing(1:6))/(E32(2)-E32(1)), E32(2:7)-E32(1:6)]
[(EDuffing(2:7) - EDuffing(1:6))*(E32(2)-E32(1)), E32(2:7)-E32(1:6)]
PhiE32 = P32\diag([0;phi32(2:end)])*P16;
[0;phi32(2:end)]
[0;phi32(2:end)]/pi
PhiE32 = P32\diag([0;phi32(2:end)])*P32;
PhiE32
PhiE32(1:7,1:7)
PhiE7 = PhiE32(1:7,1:7)
PhiE7(PhiE7<1e-9) = 0
PhiE7 = PhiE32(1:7,1:7)
PhiE7(abs(PhiE7)<1e-9) = 0
sqrt(1:7)/2
p32 = (-16:15).';
size(p32)
px32 = p2x(diag([0; px32(2:end)]))
px32 = p2x(diag([0; p32(2:end)]))
pE32 = P32'*px32*P32;
pE32
PhiE7*sqrt(5)
PhiE32(1:9,1:9)
PhiE7*sqrt(5)
px32(1:7,1:7)
-1i*px32(1:7,1:7)
-1i*pE32(1:7,1:7)
-1i*pE32(1:7,1:7)*sqrt(2/2.5)
-1i*pE32(1:9,1:9)*sqrt(2/2.5)
-1i*pE32(1:11,1:11)*sqrt(2/2.5)
E32(4:10)-E32(1:7)
pE7 = pE32(1:7,1:7)
pE7(abs(pE7)<1e-9) =
pE7(abs(pE7)<1e-9) = 0
-1i*pE7
whos
save transmonH25
E32(2:3)-E32(1:3)
E32(2:3)-E32(1:2)
(E32(2)-E32(1))-(E32(3)-E32(2))
((E32(2)-E32(1))-(E32(3)-E32(2)))/(E32(2)-E32(1))
(E32(2)-E32(1))*5/6
5/6
1/0.075
pE7
doc diag
diag(-1i*pE7,1)
diag(1i*pE7,1)
diag(1i*pE7,1)*sqrt(2/2.5)
diag(diag(1i*pE7,1)*sqrt(2/2.5), 1)
spdiags(diag(1i*pE7,1)*sqrt(2/2.5), 1, 7, 7)
spdiags([0; diag(1i*pE7,1)*sqrt(2/2.5)], 1, 7, 7)
spdiags([diag(-1i*pE7,1)*sqrt(2/2.5);0], -1, 7, 7)
spdiags([diag(1i*pE7,1)*sqrt(2/2.5);0], -1, 7, 7)
sqrt(1:7)/sqrt(2)
sqrt(1:7)
spdiags([diag(-1i*pE7,-1)*sqrt(2/2.5);0], -1, 7, 7)
aq = spdiags([0; diag(1i*pE7,1)*sqrt(2/2.5)], 1, 7, 7)
adagq = spdiags([diag(-1i*pE7,-1)*sqrt(2/2.5);0], -1, 7, 7)
save transmonH25
E7 = E32(1:7)
E7h=E7(2)*(0:6).'
Enl = E7-E7h
EnlDuffing=EDuffing-0.5-n
EnlDuff_adj = n.*(n-1)*Enl(3)/2
E7h=(E7(2)-E7(1))*(0:6).'
Enl=(E7-E7(1)) - E7h
save transmonH25