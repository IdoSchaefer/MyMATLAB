z = (0.5:0.01:1.5).';
[relE, result_direct, result_taylor, allpolydeg] = bFcomparison(z, 15);
[z, relE, result_direct, result_taylor, allpolydeg]
figure
plot(z, log10(relE))
xlabel('z')
ylabel('log(relE)')
