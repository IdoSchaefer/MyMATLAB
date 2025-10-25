function [a,w]=dvd(zint,z,w,fun,time,ro)
m=length(zint);
k=length(z);
s=length(time);
zz=[zint;z];
for ii=1:s
    for j=1:k
        w(m+j,ii)=feval(fun,z(j),time(ii));
        for i=m+j-1:-1:1
            w(i,ii)=ro*(w(i+1,ii)-w(i,ii))/(zz(m+j)-zz(i));
        end
        a(j,ii)=w(1,ii);
    end
end