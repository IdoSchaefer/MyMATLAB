function [v,eg,error,matvecs,flag]=evarnoldi(matvec,n,k,m,tol,cycles)
flag=1;
matvecs=0;
h=[];
load fr r
v=r;
error=0;
matvecs=k;
for j=1:cycles
    [v,h]=arnoldi(v,h,matvec,m);
    matvecs=matvecs+m-k;
    [q,eg]=eig(h(1:m,1:m));
    eg=diag(eg);
    [aeg,p]=sort(-abs(eg));
    z=eg(p(k+1:m));
    eg=eg(p(1:k));
    q=q(:,p(1:k));
    er=h(m+1,m)*q(m,1:k);
    error=max(abs(er));
    if error < tol
        v=v(:,1:m)*q;
        flag=0;
        return
    end
    y=qpol(h,z);
    [vs,h]=sarnoldi(y,h);
    v=v*vs;
end
