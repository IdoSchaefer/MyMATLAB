function  [z,eg,q,p]=Eleja_points(h,k,eta,ireal)
[mp1,m]=size(h);
z=[];
for j=m-3:m
    [q,eg]=eig(h(1:j,1:j));
    eg=diag(eg);
    z=[z;eg];
end
if eta==1
    [w,p]=sort(-abs(eg));
    [w1,p1]=sort(abs(z));
elseif eta==-1
     [w,p]=sort(abs(eg));
    [w1,p1]=sort(-abs(z));
else
    [w,p]=sort(abs(eg-eta));
    [w1,p1]=sort(-abs(z-eta));
end
z=z(p1);
l=length(z);
n=floor(l/2);
if ireal==1
    z=real_leja(z(1:n),k);
elseif ireal==2
    z=complex_leja(z(1:n),k);
end
eg=eg(p);
q=q(:,p);
