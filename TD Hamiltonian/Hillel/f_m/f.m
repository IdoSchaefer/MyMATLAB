function [y,k]=f(z,m)
nz=length(z);
one=ones(nz,1);
y=one;
r=y;
k=m;
while norm(r,inf)/norm(y,inf)> 1.e-13
    k=k+1;
    r=(r.*z)/k;
    y=y+r;
end