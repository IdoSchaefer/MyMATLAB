z = zeros(1, 200);
for zi=1:200
    z(zi) = rand-1 + (rand*2-1)*1i;
end
lejaz = leja(z);
plot(real(lejaz(1:50)), imag(lejaz(1:50)), '*')