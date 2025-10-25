function  [v,h]=onematvec(matvec,v,h)
[m,m1]=size(h);
v(:,m+1)=feval(matvec,v(:,m));
for i=1:m
    h(i,m)=v(:,i)'*v(:,m+1);
    v(:,m+1)=v(:,m+1)-h(i,m)*v(:,i);
end
h(m+1,m)=norm(v(:,m+1));
if h(m+1,m) < 1.e-14
    h(m+1,m)=1.e-14;
end
v(:,m+1)=v(:,m+1)/h(m+1,m);

