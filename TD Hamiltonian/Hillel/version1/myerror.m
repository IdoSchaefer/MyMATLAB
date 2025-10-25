function  [er,y,e]=myerror(h,z,a,nr,ro)
[mp1,m]=size(h);
e=zeros(mp1,1);
e(1)=nr;
[mmm,s]=size(a);
y=zeros(m,s);
flag=1;
for j=1:m
    y=y+e(1:m,1)*a(j,:);
    e(1:j+1)=(h(1:j+1,1:j)*e(1:j)-z(j)*e(1:j+1))/ro;
end
nr=norm(e);
er=max(abs(a(m,:)))*nr;

