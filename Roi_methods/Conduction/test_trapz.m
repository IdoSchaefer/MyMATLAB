minN = 5;
Nsamp = 15;
allN = zeros(1, Nsamp);
aller = zeros(1, Nsamp);
for degi = 1:Nsamp
    deg = log10(minN) + (degi-1)*0.1;
    N = round(10^deg);
    allN(degi) = N;
    h = 1/N;
    x = 0:h:1;
    f = x.^2;
    int = eq_space_trapz(f, 0, 1);
    error = int - 1/3;
    aller(degi) = error;
end
a = polyfit(log10(allN), log10(aller), 1);
display('Slope: ')
a(1)
figure
plot(log10(allN), log10(aller), 'o-')
xlabel('log(N)')
ylabel('log(error)')