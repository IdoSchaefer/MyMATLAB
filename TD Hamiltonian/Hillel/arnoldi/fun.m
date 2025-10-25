function y=fun(z,t)
global mmcheb
m=mmcheb;
zt=z*t;
if abs(zt) < m+1
    y=1;
    r=1;
    k=m;
    while abs(r)/abs(y)> 1.e-13
        k=k+1;
        r=r*zt/k;
        y=y+r;
    end
else
    y=exp(zt);
    for j=1:m
        y=j*(y-1)/zt;
     end
end
y=y*t^m;
