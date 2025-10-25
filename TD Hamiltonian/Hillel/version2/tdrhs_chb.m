function u=tdrhs(u,matvec,t,dt,m1,m2,z1,z2)
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
%          z1,z2   -   The eigenvalues of A are in the interval [z1,z2] 
%
% output:
%          u       -   matrix of dimensions  nXm1
%                       whose columns are the approximated solution at m1 Chebyshev
%                       points in the interval [t+dt,t+2*dt].
%
global c tvec tj 
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
    c=chebcoef(z1,z2,m1,m2,'fun',tvec);
    maxit=2;
    for  iter=1:maxit
        a=dvd_vec(tvec,fs);
        s=convrt(tvec,a);
        v=u(:,1);
        for j=1:m1
            v(:,j+1)=(feval(matvec,v(:,j),thalf)+s(:,j))/j;
        end
        w=chebfun(v(:,m1+1),matvec,thalf,m2,z1,z2);
        u=v(:,1:m1)*tj'+w*c;
        fs=rhs(u,t+tvec,t+dt/2);
    end
    u=new_guess(v,w,m1,m2,dt,tvec,z1,z2,'fun');
else
    thalf=t+dt/2;
    fs=rhs(u,t+tvec,thalf);
    a=dvd_vec(tvec,fs);
    s=convrt(tvec,a);
    v=u(:,1);
    for j=1:m1
        v(:,j+1)=(feval(matvec,v(:,j),thalf)+s(:,j))/j;
    end
    w=chebfun(v(:,m1+1),matvec,thalf,m2,z1,z2);
    u=new_guess(v,w,m1,m2,dt,tvec,z1,z2,'fun');
end
      


    