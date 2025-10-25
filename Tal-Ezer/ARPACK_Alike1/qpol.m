function y=qpol(h,z)
[mp1,m]=size(h);
hh=h(1:m,1:m);
y=zeros(m,1);
y(1)=1;
k=length(z);
for j=1:k
    y=hh*y-z(j)*y;
end
y=y/norm(y);    
    
