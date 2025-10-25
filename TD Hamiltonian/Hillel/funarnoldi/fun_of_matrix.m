function [u,iter,flag]=fun_of_matrix(r,matvec,fun,time,m,tol,maxit,jp)
%
% The routine approximates the vector u where  u=f(A)v.
%
%input:
%         r  -        input vector
%         matvec - routine which defines the operator ( matrix ) A
%                      it has the form -    v1=matvec(v2)  and v1,v2 should
%                      satisfy   v1=A*v2.
%         fun -      routine which defines the function f.
%         time -    vector of time points (levels) where we want to compute the solution
%         m  -       size of Krylov space
%         tol-        accuracy tolerance
%         maxit -   maximal mat-vecs allowed.
%         jp -        points in the vector v where we want to compute the solution.
%
% output:
%         u - solution vector
%         flag -     =0, solution has been computed to the desired accuracy.
%                     =1, solution has not been achieved.
%
nr=norm(r);
n=length(r);
s=length(time);
np=length(jp);
u=zeros(np,s);
iter=0;
zint=[];w=[];nr=norm(r);ro=1;
flag=1;
while iter < maxit
    iter=iter+1;
    h=[];
    [u,r,zint,nr,w,ro,flag]=funarnoldi(r,matvec,m,w,zint,fun,time,nr,u,tol,jp,ro,iter);
    if flag==0
        return
    end
end
 