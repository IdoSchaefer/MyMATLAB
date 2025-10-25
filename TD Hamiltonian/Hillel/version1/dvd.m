function [a,w]=dvd(zint,z,w,fun,time,ro)
m=length(zint);
k=length(z);
zz=[zint;z];
for j=1:k
     w(m+j,:)=feval(fun,z(j),time);
     for i=m+j-1:-1:1
         w(i,:)=ro*(w(i+1,:)-w(i,:))/(zz(m+j)-zz(i));
     end
     a(j,:)=w(1,:);
end
