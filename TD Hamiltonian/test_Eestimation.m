% Positive z:
z = (0.01:0.01:10).';
bfd = bFdirect(z, 9);
bft = bFtaylor(z, 9);
relE = abs(bfd - bft)./abs(bft);
relEestimation = eps*factorial(9)./abs(z).^9;
figure
plot(z, log10(relE), 'b')
hold on
plot(z, log10(relEestimation), 'r')

% Negative z:
zn = -(0.01:0.01:10).';
bfdn = bFdirect(zn, 9);
bftn = bFtaylor(zn, 9);
relEn = abs(bfdn - bftn)./abs(bftn);
relEestimationn = eps*factorial(9)./abs(zn).^9;
figure
plot(zn, log10(relEn), 'b')
hold on
plot(zn, log10(relEestimationn), 'r')

% Imaginary z:
zim = 1i*(0.01:0.01:10).';
bfdim = bFdirect(zim, 9);
bftim = bFtaylor(zim, 9);
relEim = abs(bfdim - bftim)./abs(bftim);
relEestimationim = eps*factorial(9)./abs(zim).^9;
figure
plot(imag(zim), log10(relEim), 'b')
hold on
plot(imag(zim), log10(relEestimationim), 'r')