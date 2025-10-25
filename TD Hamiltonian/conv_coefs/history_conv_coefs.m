Ctexp = convCtexp(9, 10);
xival = 0:0.01:1;
Mcheb = cosh((0:30).'*acosh(2*xi + 1));
Mcheb = cosh((0:30).'*acosh(2*xival + 1));
gfuns = Ctexp*Mcheb;
gfunsa = conv_gfuns_sym(9, 10, xival);
tic, gfunsa = conv_gfuns_sym(9, 10, xival);toc
tic, Ctexp = convCtexp(9, 10); Mcheb = cosh((0:30).'*acosh(2*xi + 1)); gfuns = Ctexp*Mcheb; toc
tic, Ctexp = convCtexp(9, 10); Mcheb = cosh((0:30).'*acosh(2*xival + 1)); gfuns = Ctexp*Mcheb; toc
tic, Mcheb = cosh((0:30).'*acosh(2*xival + 1)); toc
tic, for j=1:1e3, Mcheb = cosh((0:30).'*acosh(2*xival + 1)); end, toc
Mcheb_rec = cheb_pols(2*xival + 1, 30);
max(max(abs(Mcheb - Mcheb_rec)))
max(max(abs(Mcheb - Mcheb_rec)./abs(Mcheb)))
tic, for j=1:1e3, Mcheb_rec = cheb_pols(2*xival + 1, 30); end, toc
tic, for j=1:1e3, Mcheb = cosh((0:30).'*acosh(2*xival + 1)); end, toc
tic, Ctexp = convCtexp(9, 10); Mcheb = cosh((0:30).'*acosh(2*xival + 1)); gfuns = Ctexp*Mcheb; toc
tic, gfuns_temp = conv_gfuns(9, 10, xival); toc
max(max(abs(gfuns - gfuns_temp)))
tic, gfuns_temp = conv_gfuns(9, 10, xival); toc
max(max(abs(gfuns - gfuns_temp)))
(gfuns(:, end) - gfunsa(:, end))./gfunsa(:, end)
(gfuns_temp(:, end) - gfunsa(:, end))./gfunsa(:, end)
max(Mcheb)
Mcheb(:, end)
Ctexp(end, :)
Ctexp(end, :).'.*Mcheb(:, end)
Ctexp(end, :).'.*Mcheb_rec(:, end)
[Ctexp(end, :).'.*Mcheb(:, end), Ctexp(end, :).'.*Mcheb_rec(:, end)]
Ctexp(end, :).'.*Mcheb(:, end)- Ctexp(end, :).'.*Mcheb_rec(:, end)
(Ctexp(end, :).'.*Mcheb(:, end)- Ctexp(end, :).'.*Mcheb_rec(:, end))./abs(Ctexp(end, :).'.*Mcheb_rec(:, end))
[(gfuns_temp(:, end) - gfunsa(:, end))./gfunsa(:, end), (gfuns(:, end) - gfunsa(:, end))./gfunsa(:, end)]
Ctexp
Ctexp.*Mcheb(:, end).'
gfunsa
gfunsa(:, end)
max(Ctexp.*Mcheb(:, end).', [], 2)
[max(Ctexp.*Mcheb(:, end).', [], 2), gfunsa(:, end)]
C1st=convC1st(10);
figure
plot(xival, gfunsa)
min(abs(gfunsa))./gfunsa
min(abs(gfunsa(:, end)))./gfunsa(:, end)
min(abs(gfunsa), [], 2)./gfunsa(:, end)
(abs(gfunsa(:, 51)), [], 2)./gfunsa(:, 101)
abs(gfunsa(:, 51))./gfunsa(:, 101)
C1st=convC1st(10);
Mxip1 = (xi + 1).^((0:34).');
Mxip1 = (xival + 1).^((0:34).');
Mxi = xi.^((0:21).');
Mxi = xival.^((0:21).');
size(C1st)
Cfm = convCfm(14, 10);
size(Cfm)
lfuns = C1st*Mxi;
lfuns(:,end)
lfuns(2:end)./lfuns(1:end-1)
lfuns(2:end, end)./lfuns(1:end-1, end)
hfuns = Cfm*Mxip1;
hfuns(2:end, end)./hfuns(1:end-1, end)
gfuns(2:end, end)./gfuns(1:end-1, end)
gfunsa(2:end, end)./gfunsa(1:end-1, end)
abs(gfunsa(:, 51))./gfunsa(:, 101)
abs(lfuns(:, 51))./lfuns(:, 101)
abs(hfuns(:, 51))./hfuns(:, 101)
gfunsa8 = conv_gfuns_sym(8, 10, xival);
gfuns8 = conv_gfuns(8, 10, xivals);
gfuns8 = conv_gfuns(8, 10, xival);
gfuns8
gfuns8(:, end)
(gfuns8(:, end) - gfuns8a(:, end))./gfuns8a(:, end)
(gfuns8(:, end) - gfunsa8(:, end))./gfunsa8(:, end)
gfunsa8(2:end, end)./gfunsa8(1:end-1, end)
gfuns5 = conv_gfuns(5, 10, xival);
gfunsa5 = conv_gfuns_sym(5, 10, xival);
(gfuns5(:, end) - gfunsa5(:, end))./gfunsa5(:, end)
gfunsa5(2:end, end)./gfunsa5(1:end-1, end)
whos
abs(gfuns8(:, 51))./gfuns8(:, 101)
abs(gfuns5(:, 51))./gfuns5(:, 101)
abs(gfuns(:, 51))./gfuns(:, 101)
abs(hfuns(:, 51))./hfuns(:, 101)
whos
save conv_coefs
gfactors = gfunsa(2:end, end)./gfunsa(1:end-1, end);
hfactors = hfuns(2:end, end)./hfuns(1:end-1, end)
lfactors = lfuns(2:end, end)./lfuns(1:end-1, end)
save conv_coefs
load errors_harmonic
whos
est_ers99
est_ers99.texp_cheap
est_ers99.texp_exact
est_ers99.conv_cheap
est_ers99.conv_exact
gfactors(1)
[est_ers99.conv_cheap*2*gfactors(1), est_ers99.conv_exact]
[est_ers99.conv_cheap*2*gfactors(1) - est_ers99.conv_exact]
(est_ers99.conv_cheap*2*gfactors(1) - est_ers99.conv_exact)./est_ers99.conv_exact
[est_ers99.conv_cheap*2*gfactors(1), est_ers99.conv_exact, aller99]
gfactors5 = gfunsa5(2:end, end)./gfunsa5(1:end-1, end)
[est_ers55.conv_cheap*2*gfactors5(1), est_ers55.conv_exact, aller55]
gfuns7 = conv_gfuns(7, 10, xival);
gfactors7 = gfuns7(2:end, end)./gfuns7(1:end-1, end)
[est_ers79.conv_cheap*2*gfactors5(1), est_ers79.conv_exact, aller79]
[est_ers7i2.conv_cheap*2*gfactors5(2), est_ers7i2.conv_exact, aller7i2]
[est_ers7i2.conv_cheap*2*gfactors7(2), est_ers7i2.conv_exact, aller7i2]
[est_ers79.conv_cheap*2*gfactors7(1), est_ers79.conv_exact, aller79]
gfactors7
hfactors
[gfactors7 hfactors]