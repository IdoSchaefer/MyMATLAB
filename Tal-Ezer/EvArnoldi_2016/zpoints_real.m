function z=zpoints_real(h,k)
[mp1,m]=size(h);
z=[];
for j=2:m
    eg=eig(h(1:j,1:j));
    z=[z;eg];
end
[zz,p]=sort(abs(z));
z=z(p);
l=length(z);
n=floor(l/2);
z=leja_real(z(1:n,1),m-k);

