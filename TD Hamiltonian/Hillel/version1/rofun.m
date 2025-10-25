function  ro=rofun(z)
n=length(z);
zp=sum(z)/n;
ro=1;
for j=1:n
    ro=ro*abs(zp-z(j))^(1/n);
end
    
