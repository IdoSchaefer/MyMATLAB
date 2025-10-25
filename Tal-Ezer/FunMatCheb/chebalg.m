function u=chebalg(r,tpoints,matvec,fun,m,tol,parm)
%c=cheb_coef(fun,300,tpoints,tol,parm(4),parm(5));
c=cheb_coef(fun,1e4,tpoints,tol,parm(4),parm(5));
u=cheb_funmat(r,matvec,c,tol,parm(4),parm(5));