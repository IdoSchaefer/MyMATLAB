function [q,r]=myqr(w)
[n,m]=size(w);
r(1,1)=norm(w(:,1));
q(:,1)=w(:,1)/r(1,1);
for j=2:m
   q(:,j)=w(:,j);
   for i=1:j-1
      r(i,j)=q(:,i)'*q(:,j);
      q(:,j)=q(:,j)-r(i,j)*q(:,i);
   end
   r(j,j)=norm(q(:,j));
   q(:,j)=q(:,j)/r(j,j);
end