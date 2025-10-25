function [u,flag]=tdrhs(u0,matvec,t0,dt,m1,m2,tol,maxit)
% the routine solves, iteratively, problems of type
%
%    U_t = H(t)U,   t0 < t < t0+dt
%
%   input:
%           u0   -       solution vector at t=t0.
%           matvec  -  routine which performs matrix-vector multiplication
%                           where the matrix is H(t1+dt/2).
%           t0   -         initial time level.
%           dt   -         time step
%           m1  -         number of discretization points at the interval [t0  t0+dt]
%           m2  -        degree of polynomial which aproximates the residual function F.
%           tol  -         accuracy tolerance.
%           maxit -      maximal number of iterations.
%           
%  output:
%           u  -     solution vector at the  m1 discretization points at the interval [t0  t0+dt].
%           flag  -   =0  solution has been computed to the desired accuracy.
%                      =1  solution has not been computed.
%
global x
n=length(u0);
flag=0;
tvec=tpoints(dt,m1);
ub=u0;
fs=rhs0(t0+tvec,t0+dt/2);
a=dvd_vec(tvec/dt*4,fs);
for  iter0=1:maxit
%
% compute the s vectors by  Newton interpolation
%
    s=convrt(tvec/dt*4,a);
    p=1;
    for j=1:m1
        s(:,j)=p*s(:,j);
        p=p*4/dt;
    end
    if iter0==1
        ss=s;
        u=x.*u0;
        clear s
        for j=1:m1
            s(:,j)=ss(j)*(x.*u0);
        end
    end
  
%
% compute the s vectors by cosine transform.
%     a=cheb0(conj(fs'));
%     a=conj(a');
%     for i=1:n
%         s(i,:)=cheb2tay(a(i,:),dt);
%     end
    v=u0;
    tj=ones(m1,1);
    for j=1:m1
        u(:,j)=u0;
    end
    for j=1:m1-1
        v=(feval(matvec,v)+s(:,j))/j;
        tj=tj.*tvec;
        u=u+v*tj';
    end
    v=(feval(matvec,v)+s(:,m1))/m1;
    jp=1:n;
    [ww,flag]=fun_of_matrix(v,matvec,'fun',tvec,m2,tol*10,maxit,jp);
    u=u+ww;
    er=norm(u(:,m1)-ub)/norm(ub);
    [iter0 norm(u) er]
     if er < tol
         return
     end
     ub=u(:,m1);
     fs=rhs(u,t0+tvec,t0+dt/2);
     a=dvd_vec(tvec/dt*4,fs);
end


    