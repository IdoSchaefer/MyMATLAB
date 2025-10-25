function [u,tvec,flag]=tdrhs(u0,matvec,t,dt,mcheb,tol,maxit)
%
% The routine is applied for problem of type
%
%   u_t  =  Au + f(t)Bu    where A and B are square matrices.
%
%  Given u(t) as input, it computes u(t+dt).
%
% input:
%           u0 - the solution at time level t
%           matvec - routine which represents the matrix A.  matvec is of the form  
%                        u=matvec(v).  It means that u=Av.
%           t  -  the current time level.
%           dt - time interval.
%           mcheb - number of Chebyshev points in the interval [t , t+dt] and also the dimension of 
%                       Krylov space used to approximate the operator F.
%           tol - accuracy tolerance
%           maxit - limit of the number of iterations.
%
% output:
%           u - the computed solution at the interval [ t , t+dt].
%           tvec - Chebshev points at the interval [ t , t+dt ]  where the solution is computed.
%           flag  -  =0 , the solution is computed to the desired accuracy
%                     =1 , otherwise
%           
n=length(u0);
flag=0;
iu=sqrt(-1);
tvec=tpoints(dt,mcheb);
for j=1:mcheb
    u(:,j)=u0;
end
ub=u;
thalf=t+dt/2;
fs=rhs(u,t+tvec,thalf);
for  iter=1:maxit
    a=cheb(fs);
    s=cheb2tay(a,dt);
    v=u0;
    tj=ones(mcheb,1);
    for j=1:mcheb
        u(:,j)=u0;
    end
    for j=1:mcheb-1
        v=(feval(matvec,v)+s(:,j))/j;
        tj=tj.*tvec;
        u=u+v*tj';
    end
    v=(feval(matvec,v)+s(:,mcheb))/mcheb;
    jp=1:n;
    [ww,flag]=fun_of_matrix(v,matvec,'fun',tvec,mcheb,tol,maxit,jp);
    u=u+ww;
    fs=rhs(u,t+tvec,t+dt/2);
     er=norm(u-ub)/norm(ub);
     [iter norm(u) er]
     if er < tol
         return
     end
     ub=u;
end


    