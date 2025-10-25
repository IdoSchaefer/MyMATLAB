function [v,r]=gs(u)
[n,m]=size(u);
r(1,1)=norm(u(:,1));
v(:,1)=u(:,1)/r(1,1);
for j=2:m
    v(:,j)=u(:,j);
    for i=1:j-1
        r(i,j)=v(:,i)'*u(:,j);
        v(:,j)=v(:,j)-r(i,j)*v(:,i);
    end
    r(j,j)=norm(v(:,j));
    v(:,j)=v(:,j)/r(j,j);
end