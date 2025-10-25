function [u,r,zint,nr,w,ro,flag]=funarnoldi(v,matvec,m,w0,zint,fun,time,nr,u,tol,jp,ro,iter)
flag=1;
v=v/norm(v);
zp=[];
for j=1:m
    v(:,j+1)=feval(matvec,v(:,j));
    for i=1:j
        h(i,j)=v(:,i)'*v(:,j+1);
        v(:,j+1)=v(:,j+1)-h(i,j)*v(:,i);
    end
    h(j+1,j)=norm(v(:,j+1));
    v(:,j+1)=(1/h(j+1,j))*v(:,j+1);
    z=eig(h(1:j,1:j));
    zp=[zp;z];
    if j==m
        z=leja_complex(zp,zint,m,ro);
        if iter==1
            ro=rofun(z);
        end
    end
    [a,w]=dvd(zint,z,w0,fun,time,ro);
    [er,y,e]=myerror(h,z,a,nr,ro);
     [j er];
    if er < tol
        u=u+v(jp,1:j)*y;
        r=v*e;
        flag=0;
        return
    end
end
u=u+v(jp,1:j)*y;
r=v*e;
zint=[zint;z];
nr=norm(e);

