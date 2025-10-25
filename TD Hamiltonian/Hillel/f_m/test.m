function test
for j=1:10
    z=j/100;
    for m=1:10
        [y,k]=f(z,m);
        [y m k]
    end
end
j