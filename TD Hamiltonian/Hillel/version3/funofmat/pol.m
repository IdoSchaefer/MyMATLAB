function [w,r,beta]=pol(v,h,z,a,beta,ro,ireal)
[s1,s2]=size(a);
[mp1,m]=size(h);
h(:,mp1)=zeros(mp1,1);
e=zeros(mp1,1);
y=zeros(mp1,s2);
e(1)=beta;
y(1,:)=beta*a(1,:);
for j=1:m-1
    e=(1/ro)*(h*e-z(j)*e);
    y=y+e*a(j+1,:);
end
e=(1/ro)*(h*e-z(m)*e);
beta=norm(e);
e=e/beta;
if ireal==1
    w=v(:,1:m)*real(y(1:m,:));
    r=v*real(e);
else
    w=v(:,1:m)*y(1:m,:);
    r=v*e;
end

    