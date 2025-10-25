minN = 4;
Nsamp = 15;
allN = zeros(1, Nsamp);
aller = zeros(1, Nsamp);
for degi = 1:Nsamp
    deg = log10(minN) + (degi-1)*0.2;
    N = round((10^deg)/2)*2;
    allN(degi) = N;
    h = pi/N;
    x = 0:h:pi;
    f = sin(x);
    int = eq_space_simpson(f, 0, pi);
    error = int - 2;
    aller(degi) = error;
end
a = polyfit(log10(allN), log10(aller), 1);
display('Slope: ')
a(1)
figure
plot(log10(allN), log10(aller), 'o-')
xlabel('log(N)')
ylabel('log(error)')