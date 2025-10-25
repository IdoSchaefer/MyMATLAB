function u=newton(r,tpoints,matvec,fun,m,tol,restarts,ireal)
%
% the algorithm computes the vectors
%                u_i = f(A,tpoints(i))*r    where A is a matrix.
%
%  input :
%                   r =                   an input vector
%                   tpoints =       vector of points
%                   matvec =      routine which  defines implicitly the matrix A
%                   fun =              name of routine which defines the  functio f
%                   m=                    size of Krylov space.
%                   tol =             accuracy required. The algorithm stops when    ||residual||/||u|| < tol
%                   restarts= maximal number of cycles. When the desired
%                                            accuracy is achieved then the number of mat-vecs
%                                             performed is  mXcycles.
%                   ireal =               1  , if A is real
%                                           0 , if A is complex
%
%   output
%                   u  =                 the solution vectors
%

% The comments are mine (Ido).
n=length(r);
[s1,s]=size(tpoints);
if s1>s
    tpoints=tpoints';
    s=s1;
end
z=[];w=[];u=zeros(n,s); 
beta=norm(r);
r=r/beta;
% ro0 is the previous capacity; ro is the new capacity; iro is a flag, that
% is set to 0 when the value of the capacity converges. Then, it is not
% corrected any more.
ro0=1;ro=1.e13;iro=1;
m0=m;
for j=1:restarts
     m=m0;
     [v,h]=arnoldi(r,matvec,m);
     znew=newpoints(z,h,ireal);
     m=length(znew);
     % It is possible that newpoints returns more than m points, as a
     % correction for real matrices. I have to understand it more
     % precisely.
     % In such a case, it is necessary to enlarge the Hessenberg matrix and
     % the Krylov space, for the computation of u1 and r by the function
     % pol. This is the role of onematvec.
     if m > m0
         [v,h]=onematvec(matvec,v,h);
     end
     f=feval(fun,znew,tpoints);
     % The capacity ro is computed in every cycle, until its convergence.
     % Then, iro is set to 0.
     if iro==1 
         ro=capacity(znew,z);
% w is the diagonal for the computation of the new divided differences. If 
% the capacity is changed, the values of the divided differences have to be
% corrected accordingly - for every division by (z_a-z_b), it has to be 
% multiplied by ro/ro0. beta is the norm of r, which also depend on the
% capacity. Its value has to be multiplied by (ro0/r0)^(interpolation polynomial degree).
         [w,beta]=new_w_beta(w,beta,ro,ro0);
     end
     if abs(ro-ro0)/ro0 < .1
          iro=0;
     end
     ro0=ro;
     z=[z;znew]
     % a are the divided differences. w is the new diagonal.
     [a,w]=dvd(w,z,f,ro);
     % The computation of the interpolation polynomail:
     [u1,r,beta]=pol(v,h,znew,a,beta,ro,ireal);
     u=u+u1;
     er=beta*max(abs(a(m,:)))/(1+norm(u));
     [ro er]
     if er < tol
         break
    end
end

