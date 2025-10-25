function [v,h]=Earnoldi(v,h,matvec,m)
[n,k]=size(v);
for j=k:m
    v(:,j+1)=feval(matvec,v(:,j));
    for i=1:j
        h(i,j)=v(:,i)'*v(:,j+1);
        v(:,j+1)=v(:,j+1)-h(i,j)*v(:,i);
    end
    for i=1:j
        hh=v(:,i)'*v(:,j+1);
        v(:,j+1)=v(:,j+1)-hh*v(:,i);
        h(i,j)=h(i,j)+hh;
    end
    h(j+1,j)=norm(v(:,j+1));
    v(:,j+1)=(1/h(j+1,j))*v(:,j+1);
 end

