function y=fun1(z)
m=5;
nz=length(z);
one=ones(nz,1);
if norm(z,inf) < m+1
    y=one;
    r=y;
    k=m;
    while norm(r,inf)/norm(y,inf)> 1.e-13
        k=k+1;
        r=(r.*z)/k;
        y=y+r;
    end
    k-m
else
    y=exp(z);
    for j=1:m
        y=j*((y-one)./z);
    end
end
   

