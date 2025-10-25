function [unew,u0]=tdrhs(u,matvec,t,dt,m1,m2)
% 
% This routine performs one time step of a numerical algorithm for solving
%
%               u_t = Au + N(u)
%
% input:
%          u        -   in the first step, u is the initial vector.
%                    -   From the second step, u is a matrix of dimensions  nXm1
%                        whose columns are the approximated solution at m1 Chebyshev
%                        points in the interval [t,t+dt].
%          matvec -  operator which represents the matrix A.
%          t         -   the current time level.
%          dt       -   time step
%          m1      -   number of Chebyshev points in the interval [t,t+dt].
%          m2      -   Degree of Chebyshev expansion  used to approximate fun(A)v.
%
% output:
%          unew   -  the computed solution at t+dt
%          u0      -   matrix of dimensions  nXm1
%                       whose columns are the approximated solution at m1 Chebyshev
%                        points in the interval [t+dt,t+2*dt].
%
global c tvec tj A aaa bbb
[n,m]=size(u);
if m==1
%  this is first step   
    for j=2:m1
        u(:,j)=u(:,j-1);
    end
    tvec=tpoints(dt,m1);
    tj=ones(m1,1);
    for j=1:m1-1 
        tj(:,j+1)=tj(:,j).*tvec;
    end
    thalf=t+dt/2;
     fs=rhs(u,t+tvec,thalf);
    aaa=0;bbb=-sqrt(-1)*250;
    c=chebcoef(aaa,bbb,m1,m2,'fun',tvec);
    maxit=2;
    for  iter=1:maxit
        a=cheb_row(fs);
        A=sparse(cheb2tay_mat(m1,dt));
        s=a*A;
        v=u(:,1);
        for j=1:m1
            v(:,j+1)=(feval(matvec,v(:,j),thalf)+s(:,j))/j;
        end
        w=chebfun(v(:,m1+1),matvec,thalf,m2,aaa,bbb);
        u=v(:,1:m1)*tj'+w*c;
        fs=rhs(u,t+tvec,t+dt/2);
    end
    unew=u(:,m1);
    u0=new_guess(v,w,m1,m2,dt,tvec,aaa,bbb,'fun');
else
    thalf=t+dt/2;
    fs=rhs(u,t+tvec,thalf);
    a=cheb_row(fs);
    s=a*A;
    v=u(:,1);
    for j=1:m1
        v(:,j+1)=(feval(matvec,v(:,j),thalf)+s(:,j))/j;
    end
    w=chebfun(v(:,m1+1),matvec,thalf,m2,aaa,bbb);
    unew=v(:,1:m1)*tj(m1,:)'+w*c(:,m1);
    u0=new_guess(v,w,m1,m2,dt,tvec,aaa,bbb,'fun');
end
      


    