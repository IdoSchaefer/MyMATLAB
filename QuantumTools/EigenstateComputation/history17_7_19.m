load coulomb_optV240
Vf = @(x) 1 - 1./sqrt(x.^2+1)
[Efield, es_field]  = getEfield(Vabs240, xabs240, 0.1, 0.001, 1:3, [-240 240]);
whos
figure
plot(0:0.01:0.1, Efield)
plot(0:0.01:0.1, Efield(1))
plot(0:0.01:0.1, Efield(1, :))
plot(0:0.001:0.1, Efield(1, :))
hold on
plot(0:0.001:0.1, Efield(2, :))
plot(0:0.001:0.1, Efield(3, :))
figure
viewVPmiux(es_field(:, :, 1), Vabs240, xabs240, 0:0.001:0.1, x240, 0.1)
figure
x20 = -20:0.01:20;
plot(x20, Vf(x20))
plot(x20, Vf(x20)-x20*0.082)
plot(x20, Vf(x20)-x20*0.083)
plot(x20, Vf(x20)-x20*0.084)
figure
plot(x240, conj(es_field(:,83,1)).*es_field(:,83,1))
plot(x240, conj(es_field(:,82,1)).*es_field(:,82,1))
plot(x240, conj(es_field(:,81,1)).*es_field(:,81,1))
plot(x240, conj(es_field(:,82,1)).*es_field(:,82,1))
plot(x240, conj(es_field(:,83,1)).*es_field(:,83,1))
plot(x240, conj(es_field(:,84,1)).*es_field(:,84,1))
plot(x240, conj(es_field(:,85,1)).*es_field(:,85,1))
plot(x240, conj(es_field(:,86,1)).*es_field(:,86,1))
plot(x240, conj(es_field(:,90,1)).*es_field(:,90,1))
plot(x240, conj(es_field(:,101,1)).*es_field(:,101,1))
plot(x240, conj(es_field(:,82,1)).*es_field(:,82,1))
viewVPmiux(es_field(:, :, 2), Vabs240, xabs240, 0:0.001:0.1, x240, 0.1)
[Efield82, es_field82]  = getEfield(Vabs240-0.08*xabs240, xabs240, 0.01, 0.0001, 1, [-240 240]);
figure
plot(0.08:0.0001:0.09, Efield82)
fieldv= [0:0.01:0.08, 0.0801:0.001:0.09];
[Efield, es_field]  = getEfield1(Vf, miu, fieldv, n_levels, xdomain, m)
[Efield82, es_field82]  = getEfield1(Vabs240-0.08*xabs240, xabs240, fieldv, 1, [-240 240]);
figure
plot(fieldv, Efield82)
fieldv= [0:0.01:0.08, 0.0801:1e-4:0.09];
[Efield82, es_field82]  = getEfield1(Vabs240, xabs240, fieldv, 1, [-240 240]);
plot(fieldv, Efield82)
viewVPmiux(es_field82, Vabs240, xabs240, fieldv, x240, 0.1)
fieldv= [0:1e-3:0.08, 0.0801:1e-4:0.09];
[Efield82, es_field82]  = getEfield1(Vabs240, xabs240, fieldv, 1, [-240 240]);
plot(fieldv, Efield82)
viewVPmiux(es_field82, Vabs240, xabs240, fieldv, x240, 0.1)
[Efield_fine, es_field_fine]  = getEfield(Vabs240, xabs240, 0.1, 1e-4, 1:3, [-240 240]);
figure
plot(0:1e-4:0.01, Efield_fine(1,:))
plot(0:1e-4:0.1, Efield_fine(1,:))
hold on
plot(0:1e-4:0.1, Efield_fine(2,:))
plot(0:1e-4:0.1, Efield_fine(3,:))
viewVPmiux(es_field_fine, Vabs240, xabs240, 0:1e-4:0.1, x240, 0.01)
viewVPmiux(es_field_fine(:, :, 1), Vabs240, xabs240, 0:1e-4:0.1, x240, 0.01)
viewVPmiux(es_field_fine(:, :, 2), Vabs240, xabs240, 0:1e-4:0.1, x240, 0.01)
figure
plot(x20, Vf(x20)-x20*0.03)
whos
Nx
[Efield40, es_field40]  = getEfield(Vabs, xabs, 0.1, 0.001, 1:3, xdomain);
xdomain
figure
plot(0:0.001:0.1, Efield40(1, :))
hold on
plot(0:0.001:0.1, Efield40(2, :))
plot(0:0.001:0.1, Efield40(3, :))
plot(x20, Vf(x20)-x20*0.009)
plot(x20, Vf(x20)-x20*0.036)
viewVPmiux(es_field40(:, :, 1), Vabs, xabs, 0:1e-3:0.1, x, 0.01)
viewVPmiux(es_field40(:, :, 1), Vabs, xabs, 0:1e-3:0.1, x, 0.1)
viewVPmiux(es_field40(:, :, 2), Vabs, xabs, 0:1e-3:0.1, x, 0.1)
viewVPmiux(es_field40(:, :, 3), Vabs, xabs, 0:1e-3:0.1, x, 0.1)
figure
plot(0:0.036, Efield40(2) - Efield40(1))
plot(0:0.036, Efield40(2,:) - Efield40(1,:))
plot(0:0.01:0.036, Efield40(2,1:37) - Efield40(1,1:37))
plot(0:0.001:0.036, Efield40(2,1:37) - Efield40(1,1:37))
plot(x20, Vf(x20)-x20*0.1)
plot(x20, Vf(x20)-x20*0.12)
plot(x20, Vf(x20)-x20*0.15)
[Efield40, es_field40]  = getEfield(Vabs, xabs, 0.15, 0.001, 1:3, xdomain);
figure
plot(0:1e-3:0.15, Efield40(1, :))
plot(x20, Vf(x20)-x20*0.17)
[Efield40, es_field40]  = getEfield(Vabs, xabs, 0.17, 0.001, 1:3, xdomain);
plot(0:1e-3:0.17, Efield40(1, :))
viewVPmiux(es_field40(:, :, 1), Vabs, xabs, 0:1e-3:0.17, x, 0.1)
whos
Kx = p2x(diag(K));
H01 = diag(Vabs) + Kx;
H01 = diag(Vabs-0.1*xabs) + Kx;
[P01, D01] = eig(H01);
E01 = diag(D01);
[E01, orderE01] = sort(real(E01));
P01 = P01(:, orderE01);
figure
plot(0:255, E01*ones(1, 256))
plot(0:255, E01)
plot(0:255, E01, 'o')
figure
plot(x, conj(P01(:,50)).*P01(:,50))
plot(x, conj(P01(:,49)).*P01(:,49))
plot(x, conj(P01(:,51)).*P01(:,51))
plot(x, conj(P01(:,52)).*P01(:,52))
pl
plot(x, conj(P01(:,52)).*P01(:,52))
plot(x, conj(P01(:,53)).*P01(:,53))
plot(x, conj(P01(:,54)).*P01(:,54))
plot(x20, Vf(x20)-x20*0.1)
whos
[~, ~, ~, ~, P240, ~] = gsV(Vf, xdomain240, Nx240);
XE240 = P'*diag(x240)*P240;
XE240 = P240'*diag(xabs240)*P240;
XE240(1:5, 1:5)
miu240 = [0 XE240(1,2); XE240(2,1) 0]
miu240 = real([0 XE240(1,2); XE240(2,1) 0])
H0 = diag(E240(1:2))
H0 = real(diag(E240(1:2)))
[EfieldTLS, es_fieldTLS]  = getEfieldH(H0, miu240, 0:0.01:0.2, 1:2);
figure
[EfieldTLS, es_fieldTLS]  = getEfieldH(H0, miu240, 0:0.001:0.2, 1:2);
plot(0:1e-3:0.2, EfieldTLS)
hold on
plot(0:1e-3:0.17, Efield40(1, :))
miu240_3 = real(XE240(1:3,1:3))
miu240_3(abs(miu240_3)<1e-13) = 0
H03 = real(diag(E240(1:3)))
[Efield3LS, es_field3LS]  = getEfieldH(H03, miu240_3, 0:0.001:0.2, 1:3);
figure
plot(0:1e-3:0.2, Efield3LS)
hold on
plot(0:1e-3:0.17, Efield40(1, :))
plot(0:1e-3:0.17, Efield40(2, :))
plot(0:1e-3:0.036, Efield40(2, 1:37))
[psi3, mniter, matvecs, max_errors] = SemiGlobalH(@(psi) H03*psi - miu240_3*psi*fieldfun(t).', @(psi1, t1, psi2, t2) -miu240_3*psi1*(fieldfun(t1).' - fieldfun(t2)), 1, [], [], [1;0;0], 0:0.2:1e3, 5e3, 7, 7, 1e-7);
t = 0:0.2:1e3;
w = 0:pi/1e3:pi/0.2;
dw = pi/1e3;
dctfactor = 1e3/(sqrt(5e3*pi));
load coulomb_optV240
cw = 0.1*sin(0.06*(t-500));
cww = dctI(cw)*dctfactor;
cwwfil = cww.*exp(-(w-0.06).^2/(2*0.01^2));
cwfil = dctI(cwwfil)/dctfactor;
cww_con = fieldw20b(cwwfil.', 5e5*exp(-(w.'-0.06).^2/(2*0.01^2)), dw).';
cw_con_allt = dctIintgrid(cww_con, 1e3, t_ts(1:6))/dctfactor;
Nt_ts = 7;
tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
t_ts = 0.5*(tcheb+1)*dt;
dt=0.2;
tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
t_ts = 0.5*(tcheb+1)*dt;
t_ts
cw_con_allt = dctIintgrid(cww_con, 1e3, t_ts(1:6))/dctfactor;
cw_con_window_allt = cw_con_allt*0.5.*(tanh((allt-120)/30) - tanh((allt-880)/30));
allt = [kron(0:0.2:999.8, ones(1,6)) + kron(ones(1, 5000), t_ts(1:6)), T]'
allt = [kron(0:0.2:999.8, ones(1,6)) + kron(ones(1, 5000), t_ts(1:6)), 1e3];
allt(1:30)
size(allt)
cw_con_window_allt = cw_con_allt*0.5.*(tanh((allt-120)/30) - tanh((allt-880)/30));
[allpsi3, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H03, [0 1], miu40_3, [1;0;0], [0 T], 5e3, 7, 7, 1e-5, cw_con_window_allt);
[allpsi3, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H03, [0 1], miu240_3, [1;0;0], [0 t], 5e3, 7, 7, 1e-5, cw_con_window_allt);
figure
plot(allt, cw_con_window_allt)
[allpsi3, field1, mniterc] = solveOCMLSih(@ihfieldMLS, H03, [0 1], miu240_3, [1;0;0], [0 1e3], 5e3, 7, 7, 1e-5, cw_con_window_allt);
figure
plot(allt, allpsi3.*conj(allpsi3))
plot(allt, allpsi3(:,1).*conj(allpsi3(:,1)))
plot(allt, allpsi3(1,:).*conj(allpsi3(1,:)))
size(allpsi3)
plot(t, allpsi3(1,1,:).*conj(allpsi3(1,1,:)))
psi3(:) = allpsi3(1,1,:).*conj(allpsi3(1,1,:));
size(psi3)
psi3(:,:) = allpsi3(1,:,:).*conj(allpsi3(1,:,:));
psi3(:,:) = allpsi3(:,1,:).*conj(allpsi3(:,1,:));
clear psi3
psi3(:,:) = allpsi3(:,1,:).*conj(allpsi3(:,1,:));
size(psi3)
psi3 = [psi3 allpsi3(:, 7, end)];
size(psi3)
plot(t, psi3.*conj(psi3))
npai3 = sqnorm(psi3);
figure
plot(t, npsi3)
plot(t, pai3)
plot(t, npai3)
norm(psi3(:,end))
clear npai3
edit sqnorm
edit normU
npsi3 = sqnorm(psi3);
size(npsi3)
plot(t, npsi3)
psi3(:,:) = allpsi3(:,1,:);
clear psi3;
psi3(:,:) = allpsi3(:,1,:);
psi3 = [psi3 allpsi3(:, 7, end)];
plot(t, psi3.*conj(psi3))
hold on
plot(allt, cw_con_window_allt)
whos
doc fmin
threshold_eps
1e2/(1 + 1e2^2)^(3/2) - field
1e-2/(1 + 1e-2^2)^(3/2) - field
1e2/(1 + 1e2^2)^(3/2) - 1e-3
0.5/(1 + 0.5^2)^(3/2) - 1e-3
threshold_eps
figure
plot(1e-3:1e-3:0.1, threshold_field)
threshold_eps
plot(1e-3:1e-3:0.1, threshold_field)
figure
plot(1e-3:1e-3:0.1, xthreshold_field)
hold on
plot(1e-3:1e-3:0.1, Efield40(:, 2:101))
figure
plot(1e-3:1e-3:0.1, threshold_field - Efield40(1, 2:101))
plot(1e-3:1e-3:0.036, threshold_field - Efield40(1, 2:37))
plot(1e-3:1e-3:0.036, threshold_field(1:37) - Efield40(1, 2:37))
plot(1e-3:1e-3:0.036, threshold_field(2:37) - Efield40(1, 2:37))
plot(1e-3:1e-3:0.1, threshold_field - Efield40(1, 2:101))
hold on
plot(1e-3:1e-3:0.036, threshold_field(2:37) - Efield40(1, 2:37))
plot(1e-3:1e-3:0.036, Efield40(2, 2:37) - Efield40(1, 2:37))
plot(0:1e-3:0.036, Efield40(2, 1:37) - Efield40(1, 1:37))
gs_oc_endw3
gs_oc_endw_sec3
1 -
1-allgsocsecw
load('gsoc_comparison.mat')
1-allgsocsecw
1-allgsocsecw3
figure
plot(0:0.2:1e3, allfields_psi3_sec(:,:,1))
plot(0:0.2:1e3, allfields_psi3_sec(:,:,1).*conj(allfields_psi3_sec(:,:,1)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,1).*conj(allfields_psi3_sec(1,:,1)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,2).*conj(allfields_psi3_sec(1,:,2)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,3).*conj(allfields_psi3_sec(1,:,3)))
figure
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,1).*conj(allfields_psi3_sec(1,:,1)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,2).*conj(allfields_psi3_sec(1,:,2)))
figure
plot(allt, cw_con_window_allt_sec)
figure
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,4).*conj(allfields_psi3_sec(1,:,4)))
figure
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,4).*conj(allfields_psi3_sec(1,:,5)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,5).*conj(allfields_psi3_sec(1,:,5)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,6).*conj(allfields_psi3_sec(1,:,6)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,7).*conj(allfields_psi3_sec(1,:,7)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,8).*conj(allfields_psi3_sec(1,:,8)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,9).*conj(allfields_psi3_sec(1,:,9)))
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,10).*conj(allfields_psi3_sec(1,:,10)))
plot(0:0.2:1e3, allfields_psi3_sec(:,:,5).*conj(allfields_psi3_sec(:,:,5)))
plot(0:0.2:1e3, allfields_psi3_sec(2:3,:,10).*conj(allfields_psi3_sec(2:3,:,10)))
figure
plot(0:0.2:1e3, 1-allfields_psi3(:,:,5).*conj(allfields_psi3(:,:,5)))
plot(0:0.2:1e3, 1-allfields_psi3(1,:,5).*conj(allfields_psi3(1,:,5)))
[Efield3LS005, es_field3LS005]  = getEfieldH(H03, miu240_3, cw_con_window_allt*0.5, 1:3);
gsoc005 = es_field3LS005(:,:,1)'*allfields_psi3(:, :, 5);
gsoc005 = diag(gsoc005);
size(gsoc005)
hold on
gsoc005 = gsoc005.*conj(gsoc005);
plot(0:0.2:1e3, gsoc005)
plot(0:0.2:1e3, 1-gsoc005)
figure
plot(0:0.2:1e3, Efield005(0:0.2:1e3))
plot(0:0.2:1e3, Efield005)
plot(0:0.2:1e3, Efield3LS005)
plot(0:0.2:1e3, Efield3LS005(1,:))
size(Efield3LS005)
size(es_field3LS005)
size(allfields_psi3)
gsoc005 = es_field3LS005(:,1:6:end,1)'*allfields_psi3(:, :, 5);
size(gsoc005)
gsoc005 = diag(gsoc005);
gsoc005 = gsoc005.*conj(gsoc005);
plot(0:0.2:1e3, 1-gsoc005)
figure
plot(0:0.2:1e3, (1-gsoc005)./1-allfields_psi3(1,:,5).*conj(allfields_psi3(1,:,5)))
plot(0:0.2:1e3, (1-gsoc005)./(1-allfields_psi3(1,:,5).*conj(allfields_psi3(1,:,5))))
size(gsoc005)
plot(0:0.2:1e3, (1-gsoc005.')./(1-allfields_psi3(1,:,5).*conj(allfields_psi3(1,:,5))))
[Efield3LS005, es_field3LS005]  = getEfieldH(H03, miu240_3, cw_con_window_allt*0.5, 1:3);
[Efield3LS, es_field3LS]  = getEfieldH(H03, miu240_3, 0:0.001:0.2, 1:3);
figure
plot(0:1e-3:0.2, Efield3LS)
gsoc005 = es_field3LS005(:,1:6:end,1)'*allfields_psi3(:, :, 5);
gsoc005 = diag(gsoc005);
gsoc005 = gsoc005.*conj(gsoc005);
plot(0:0.2:1e3, 1-gsoc005)
[Efield3LS01, es_field3LS01]  = getEfieldH(H03, miu240_3, cw_con_window_allt, 1:3);
gsoc01 = es_field3LS01(:,1:6:end,1)'*allfields_psi3(:, :, 5);
gsoc01 = diag(gsoc01);
gsoc01 = gsoc01.*conj(gsoc01);
figure
plot(0:0.2:1e3, 1-allfields_psi3(1,:,10).*conj(allfields_psi3(1,:,10)))
hold on
plot(0:0.2:1e3, 1-gsoc01)
gsoc01 = es_field3LS01(:,1:6:end,1)'*allfields_psi3(:, :, 10);
gsoc01 = diag(gsoc01);
gsoc01 = gsoc01.*conj(gsoc01);
plot(0:0.2:1e3, 1-gsoc01)
4.55e-6*93143.7653
4.55e-6*159855.9745
4.55e-6*134041.8400
4.55e-6*93143.7653
4.55e-6*79971.7422
4.55e-6*67067.547
45.6/780]
45.6/780
fundamental_rare_gases = [7.2734e-01, 6.0989e-01, 4.2380e-01, 3.6387e-01
, 3.0516e-01]
fundamental_rare_gases = [7.2734e-01, 6.0989e-01, 4.2380e-01, 3.6387e-01, 3.0516e-01]
TiSapphire_w = 5.8462e-02
fundamental_rare_gases/TiSapphire_w
fundamental_rare_gases/(2*TiSapphire_w)
save fundamental_rare_gases fundamental_rare_gases TiSapphire_w
fundamental_rare_gases/(3*TiSapphire_w)
hold on
plot(allt, cw_con_window_allt)
viewVPmiux(es_field40(:, :, 1), Vabs, xabs, 0:1e-3:0.17, x, 0.1)
set(Pcurve, 'Ydata', P(:, 101)), set(Vcurve, 'Ydata', V(:, 101))
figure
plot(x, es_field40(:, 101, 1).*conj(es_field40(:, 101, 1)))
plot(x, es_field40(:, 51, 1).*conj(es_field40(:, 51, 1)))
0.232/0.06
0.395/0.12
1-allgsoc
(0.07/5.33e-9)^2
(0.06/5.33e-9)^2
(0.05/5.33e-9)^2
plot(x, es_field40(:, 61, 1).*conj(es_field40(:, 61, 1)))
plot(x, es_field40(:, 71, 1).*conj(es_field40(:, 71, 1)))
[Efield3LSsec003, es_field3LSsec003]  = getEfieldH(H03, miu240_3, cw_con_window_allt_sec*0.3, 1:3);
1-allgsocsecw3
figure
plot(0:0.2:1e3, Efield3LSsec003(1,:))
size(Efield3LSsec003)
plot(allt, Efield3LSsec003)
figure
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,3).*conj(allfields_psi3_sec(1,:,3)))
gsocsec003 = es_field3LSsec003(:,1:6:end,1)'*allfields_psi3_sec(:, :, 3);
gsocsec003 = diag(gsocsec003);
gsocsec003 = gsocsec003.*conj(gsocsec003);
hold on
plot(0:0.2:1e3, 1-gsocsec003)
figure
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,4).*conj(allfields_psi3_sec(1,:,4)))
gsocsec004 = es_field3LSsec003(:,1:6:end,1)'*allfields_psi3_sec(:, :, 4);
gsocsec004 = es_field3LSsec004(:,1:6:end,1)'*allfields_psi3_sec(:, :, 4);
[Efield3LSsec004, es_field3LSsec004]  = getEfieldH(H03, miu240_3, cw_con_window_allt_sec*0.4, 1:3);
gsocsec004 = es_field3LSsec004(:,1:6:end,1)'*allfields_psi3_sec(:, :, 4);
gsocsec004 = gsocsec004.*conj(gsocsec004);
gsocsec004 = es_field3LSsec004(:,1:6:end,1)'*allfields_psi3_sec(:, :, 4);
gsocsec004 = diag(gsocsec004);
gsocsec004 = gsocsec004.*conj(gsocsec004);
figure
hold on
plot(0:0.2:1e3, 1-gsocsec004)
1-allgsocsec
[~, ~, psisec003] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), 0.1*3*cww_con_window, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
figure
ocsec003 = P240'*psisec003;
ocsec003 = ocsec003.*conj(ocsec003);
figure
plot(0:0.2:1e3, ocsec003)
plot(0:0.2:1e3, 1-ocsec003(1,:))
[~, ~, psisec003] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), 0.1*3*cww_consec_window, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
ocsec003 = P240'*psisec003;
ocsec003 = ocsec003.*conj(ocsec003);
plot(0:0.2:1e3, 1-ocsec003(1,:))
[Efield40, es_field40]  = getEfield(Vabs, xabs, 0.17, 0.001, 1:3, xdomain);
[Efield40sec003, es_field40sec003]  = getEfield1(Vabs, xabs, cw_consec_window, 1:10, xdomain);
x(32)
x(33)
x(65)
x(129)
x(193)
gsocsec004 = es_field40sec003(65:193,:,1)'*psisec003(:, :);
x240(392)
x240(384)
x240(385)
385-65
x240(320)
x240(321)
385+64
x240(449)
gsocsec004 = es_field40sec003(65:193,:,1)'*psisec003(:, 385:449);
size(es_field40sec003)
size(psisec003)
gsocsec004 = es_field40sec003(65:193,:,1)'*psisec003(:, 321:449);
gsocsec004 = es_field40sec003(65:193,:,1)'*psisec003(321:449,:);
figure
hold on
plot(0:0.2:1e3, 1-gsocsec004)
clf
gsocsec004 = diag(gsocsec004);
gsocsec004 = gsocsec004.*conj(gsocsec004);
plot(0:0.2:1e3, 1-allfields_psi3_sec(1,:,4).*conj(allfields_psi3_sec(1,:,4)))
plot(0:0.2:1e3, 1-ocsec003(1,:))
plot(0:0.2:1e3, 1-gsocsec004)
gsocsec003 = gsocsec004;
clear gsocsec004
plot(0:0.2:1e3, 1-ocsec003(1,:))
hold on
plot(0:0.2:1e3, 1-gsocsec003)
1-allgsocsecw
figure
[Efield40sec003, es_field40sec003]  = getEfield1(Vabs, xabs, cw_consec_window(1:5:end)*0.3, 1:10, xdomain);
gsocsec003 = es_field40sec003(65:193,:,1)'*psisec003(321:449,:);
gsocsec003 = diag(gsocsec003);
gsocsec003 = gsocsec003.*conj(gsocsec003);
plot(0:0.2:1e3, 1-gsocsec003)
gsocsec003 = es_field40sec003(65:193,:,1)'*psisec003(321:449,1:5:end);
gsocsec003 = diag(gsocsec003);
gsocsec003 = gsocsec003.*conj(gsocsec003);
size(gsocsec003)
plot(0:1:1e3, 1-gsocsec003)
figure
plot(x240, P240(:,1).*conj(P240(:,1)))
plot(x240, P240(:,2).*conj(P240(:,2)))
plot(x240, P240(:,3).*conj(P240(:,3)))
miu240_3
P240(:,2)'*(x240.*P240(:, 1))
P240(:,3)'*(x240.*P240(:, 2))
1-allgsocsecw3
[~, ~, psisec004] = guessresults_pnaE0b(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.5*(tanh(50*(x-0.9)) - tanh(5)), 0.1*4*cww_consec_window, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), 0, 1e3, 0.2, 7, 7, 1e-4);
[Efield40sec004, es_field40sec004]  = getEfield1(Vabs, xabs, cw_consec_window(1:5:end)*0.3, 1:3, xdomain);
plot(1e-3:0.009, Efield40(3, 1:9))
Efield40(3, 1:9)
plot(1e-3:0.009, real(Efield40(3, 1:9)))
plot(1e-3:1e-3:0.009, real(Efield40(3, 1:9)))
plot(0:1e-3:1e-3:0.009, real(Efield40(3, 1:10)))
plot(0:1e-3:0.009, real(Efield40(3, 1:10)))
plot(0:1e-3:0.009, Efield40(3, 1:10) - Efield40(1, 1:10))
[Efield40sec004, es_field40sec004]  = getEfield1(Vabs, xabs, cw_consec_window(1:5:end)*0.4, 1:3, xdomain);
gsocsec004 = es_field40sec004(65:193,:,1)'*psisec004(321:449,1:5:end);
gsocsec004 = diag(gsocsec004);
gsocsec004 = gsocsec004.*conj(gsocsec004);
figure
ocsec004 = P240'*psisec004;
ocsec004 = ocsec004.*conj(ocsec004);
figure
plot(0:0.2:1e3, 1-ocsec004(1,:))
1-allgsocsecw
hold on
plot(0:1:1e3, 1-gsocsec004)
plot(0:0.2:1e3, 1-ocsec004(2,:))
plot(0:0.2:1e3, ocsec004(2,:))
plot(0:0.2:1e3, ocsec004(3,:))
size(ocsec004)
plot(0:0.2:1e3, ocsec004(4,:))
plot(0:0.2:1e3, ocsec004(5,:))
figure
plot(E240, ocsec004(:, end), '*')
plot(E240(2:end), ocsec004(2:end, end), '*')
figure
plot(2:768, ocsec004(2:end, end), '*')
figure
plot(E240(36:90), ocsec004(36:90))
plot(E240(36:90), ocsec004(36:90, end))
E240(1:10)
figure
plot(E240(91:191), ocsec004(91:191, end))
figure
plot(E240(2:150), ocsec004(2:150, end))
plot(0:0.2:1e3, ocsec004(4,:))
plot(0:0.2:1e3, ocsec004(42,:))
plot(0:0.2:1e3, sum(ocsec004(36:48,:)))
plot(0:0.2:1e3, sum(ocsec004(36:2:48,:)))
plot(0:0.2:1e3, sum(ocsec004(83:2:91,:)))
figure
plot(0:0.2:1e3, allfields_psi3_sec(2:3,:,4).*conj(allfields_psi3_sec(2:3,:,4)))
figure
plot(x240, P240(:,42).*conj(P240(:,42)))
plot(x240, P240(:,4).*conj(P240(:,4)))
hold on
plot(x240, P240(:,42).*conj(P240(:,42)))
plot(x240, P240(:,36).*conj(P240(:,36)))
plot(x240, P240(:,38).*conj(P240(:,38)))
plot(x240, P240(:,40).*conj(P240(:,40)))
plot(x240, P240(:,42).*conj(P240(:,42)))
plot(x240, P240(:,44).*conj(P240(:,44)))
figure
plot(0:0.2:1e3, allfields_psi3_sec(2:3,:,3).*conj(allfields_psi3_sec(2:3,:,3)))
plot(0:1:1e3, 1-gsocsec003)
plot(0:0.2:1e3, allfields_psi3_sec(2:3,:,3).*conj(allfields_psi3_sec(2:3,:,3)))
hold on
plot(0:1:1e3, 1-gsocsec003)
whos
gsocsec3LS003 = es_field3LSsec00(:,:,3)'*allfields_psi3_sec(:, :, 3);
gsocsec3LS003 = es_field3LSsec003(:,:,3)'*allfields_psi3_sec(:, :, 3);
size(gsocsec3LS003)
gsocsec3LS003 = es_field3LSsec003(:,1:6:end,3)'*allfields_psi3_sec(:, :, 3);
size(gsocsec3LS003)
gsocsec3LS003 = diag(gsocsec3LS003);
gsocsec3LS003 = gsocsec3LS003.*conj(gsocsec3LS003);
plot(0:0.2:1e3, 1-gsocsec3LS003)
gsocsec3LS003 = es_field3LSsec003(:,1:6:end,1)'*allfields_psi3_sec(:, :, 3);
gsocsec3LS003 = diag(gsocsec3LS003);
gsocsec3LS003 = gsocsec3LS003.*conj(gsocsec3LS003);
plot(0:0.2:1e3, 1-gsocsec3LS003)
gsocsec3LS004 = es_field3LSsec004(:,1:6:end,1)'*allfields_psi3_sec(:, :, 4);
gsocsec3LS004 = diag(gsocsec3LS004);
gsocsec3LS004 = gsocsec3LS004.*conj(gsocsec3LS004);
plot(0:0.2:1e3, 1-gsocsec3LS004)
plot(0:0.2:1e3, allfields_psi3_sec(2:3,:,3).*conj(allfields_psi3_sec(2:3,:,3)))
plot(0:0.2:1e3, allfields_psi3_sec(2:3,:,4).*conj(allfields_psi3_sec(2:3,:,4)))
hold on
plot(0:0.2:1e3, 1-gsocsec3LS004)