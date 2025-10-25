function u=myfun(d,p,invp,fun,v)
dd=feval(fun,d);
u=p*(dd.*(invp*v));