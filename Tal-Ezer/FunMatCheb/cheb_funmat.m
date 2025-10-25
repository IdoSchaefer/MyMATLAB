function u=cheb_funmat(v,matvec,c,tol,a,b)
m=length(c);
p1=b-a;
p2=b+a;
t1=v; 
t2=1/(b-a)*(2*feval(matvec,v)-(a+b)*v);
u=c(1)*t1+c(2)*t2;
for k=3:m
    t3=2/(b-a)*(2*feval(matvec,t2)-(a+b)*t2)-t1;
    er=norm(c(k)*t3)/(1+norm(u));
    u=u+c(k)*t3;
    if er < tol
        return
    end
    t1=t2;
    t2=t3;
end

