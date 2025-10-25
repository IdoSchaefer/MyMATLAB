function u=newton(r,tpoints,matvec,fun,m,tol,restarts,parm)
%
% the algorithm computes the vectors
%                u_i = f(A,tpoints(i))*r    where A is a matrix.
%
%  input :
%                   r =               an input vector
%                   tpoints =       vector of points
%                   matvec =      routine which  defines implicitly the matrix A
%                   fun =            name of routine which defines the  function f
%                   m=               size of Krylov space.
%                   tol =             accuracy required. The algorithm stops when    ||residual||/||u|| < tol
%                   restarts=       maximal number of cycles. When the desired
%                                      accuracy is achieved then the number of mat-vecs
%                                      performed is  m X restarts.
%                   parm(1) =     1  , if A is real
%                                     0 ,  if A is complex
%                   parm(2) =     method for computing the Ritz values
%                                     =1 eigenvalues of the Hessenberg  matrix
%                                     =2 square root of the eigenvalues of the square of the
%                                     modified Hessenberg matrix. Usually parm(2) should be 1.
%                    parm(3)      =1 ,  use Newton algorithm.
%                                    =2 , use Chebyshev algorithm.
%                                     
%
%   output
%                   u  =           the solution vectors
%
n=length(r);
s=length(tpoints);
z=[];w=[];u=zeros(n,s); 
beta=norm(r);
r=r/beta;
ro0=1;ro=1.e13;iro=1;
m0=m;
for j=1:restarts
     m=m0;
     [v,h]=arnoldi(r,matvec,m);
     znew=newpoints(z,h,parm(1),parm(2));
     m=length(znew);
     if m > m0
         [v,h]=onematvec(matvec,v,h);
     end
     clear f
     for jj=1:length(tpoints)
         f(:,jj)=feval(fun,tpoints(jj)*znew);
     end
     if iro==1 
         ro=capacity(znew,z);
         [w,beta]=new_w_beta(w,beta,ro,ro0);
     end
     if abs(ro-ro0)/ro0 < .1
          iro=0;
     end
     ro0=ro;
     z=[z;znew];
     [a,w]=dvd(w,z,f,ro);
     [u1,r,beta]=polf(v,h,znew,a,beta,ro,parm(1));
     u=u+u1;
     er=beta*max(abs(a(m,:)))/(1+norm(u));
     if er < tol
        break
    end
end

