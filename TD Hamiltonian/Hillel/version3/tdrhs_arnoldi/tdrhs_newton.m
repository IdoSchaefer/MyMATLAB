function u=tdrhs_newton(u,matvec,t,dt,m1,m2,ireal)
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
%          ireal    -    =1 the computations should be real
%                         =0 the computations  can be complex
% output:
%          u       -   matrix of dimensions  nXm1
%                       whose columns are the approximated solution at m1 Chebyshev
%                       points in the interval [t+dt,t+2*dt].
%
global tvec 
global mglobal
global un
mglobal=m1;
[n,m]=size(u);
tol=1.e-5;cycles=20;
tvec=tpoints(dt,m1);
if m==1
    for j=2:m1
        u(:,j)=u(:,j-1);
    end
    un=u(:,1);
    tj=ones(m1,1);
    for j=1:m1-1 
        tj(:,j+1)=tj(:,j).*tvec;
    end
    fs=rhs(u,un,tvec);
    maxit=4;
    for  iter=1:maxit
        un=u(:,(m1+1)/2);
        a=dvd_vec(tvec,fs);
        s=convrt(tvec,a);
        v=u(:,1);
        u=advance(v,matvec,m1,m2,tvec,s,tol,cycles,ireal);
        fs=rhs(u,un,tvec);
    end
    u=advance(v,matvec,m1,m2,tvec+dt,s,tol,cycles,ireal);
else
    un=u(:,(m1+1)/2);
    fs=rhs(u,un,t+tvec);
    a=dvd_vec(tvec,fs);
    s=convrt(tvec,a);
    v=u(:,1);
    u=advance(v,matvec,m1,m2,tvec+dt,s,tol,cycles,ireal);
end
      


    