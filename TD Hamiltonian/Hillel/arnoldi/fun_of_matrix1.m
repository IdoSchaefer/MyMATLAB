function [u,iter,flag]=fun_of_matrix1(v,matvec,fun,time,m,tol,maxit,jp)
%
% The routine approximates the vector u where  u=f(A)v.
%
%input:
%         v  -        input vector
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
%         iter -      number of mat-vecs performed
%         flag -     =0, solution has been computed to the desired accuracy.
%                     =1, solution has not been achieved.
%
nb=norm(v);
nr=norm(v);
r=v;
n=length(v);
s=length(time);
np=length(jp);
u=zeros(np,s);
iter=0;
ro = 0;
zint=[];w=[];
flg=1;
while iter < maxit
    [v,h]=arnoldi(r,matvec,m);
     [z, ro]=leja_points1(h,zint,ro);
  %  z=eig(h(1:m,1:m));
%      z=sqrt(-1)*imag(z);
%     if iter==0
%         ro=rofun(z);
%     end
    [a,w]=dvd(zint,z,w,fun,time,ro);
    zint=[zint;z];
    [u,r,nr,flag,er]=gresidum(u,v,h,z,a,nr,tol,ro,jp);
    if flag==0 
%           plot(real(zint),imag(zint),'k+')
%           axis equal
          return
    end
    iter=iter+m;
end
%  plot(real(zint),imag(zint),'k+')
%  axis equal
 