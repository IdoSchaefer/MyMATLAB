% The results are identical for a real and positive z:
z1 = (1:100).';
[relE1, result_direct1, result_taylor1, allpolydeg1] = bFcomparison(z1, 9);
[z1, relE1, result_direct1, result_taylor1, allpolydeg1]
figure
plot(z1, log10(relE1))
xlabel('z1')
ylabel('log(relE1)')
% The Taylor expansion does not converge monotonically for negative z, and
% becomes unstable for large negative z:
z2 = -(1:100).';
[relE2, result_direct2, result_taylor2, allpolydeg2] = bFcomparison(z2, 9);
[z2, relE2, result_direct2, result_taylor2, allpolydeg2]
figure
plot(z2, log10(relE2))
xlabel('z2')
ylabel('log(relE2)')
% It is evident from the following figure that the problem is with the Taylor form:
figure
plot(z2(1:64), result_direct2(1:64), 'b')
hold on
plot(z2(1:64), result_taylor2(1:64), 'r')
xlabel('z2')
ylabel('result2')
% The Taylor expansion does not converge monotonically for an imaginary z, and
% becomes unstable for large z:
z3 = 1i*(1:100).';
[relE3, result_direct3, result_taylor3, allpolydeg3] = bFcomparison(z3, 9);
[z3, relE3, result_direct3, result_taylor3, allpolydeg3]
figure
plot(imag(z3), log10(relE3))
xlabel('Im(z3)')
ylabel('log(relE3)')
% It is evident from the following figure that the problem is with the Taylor form:
figure
plot(imag(z3(1:64)), real(result_direct3(1:64)), 'b')
hold on
plot(imag(z3(1:64)), real(result_taylor3(1:64)), 'r')
xlabel('Im(z3)')
ylabel('Re(result3)')