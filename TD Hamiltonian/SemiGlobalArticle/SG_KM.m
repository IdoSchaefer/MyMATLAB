load article_data579
plot(log10(allmv5), log10(aller5), '-o')
hold on
plot(log10(allmv7), log10(aller7), '-o')
plot(log10(allmv9), log10(aller9), '-o')
xlabel('log(matvecs)')
ylabel('log(error)')
