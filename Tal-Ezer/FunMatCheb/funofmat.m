function u=funofmat(r,tpoints,matvec,fun,m,tol,restarts,parm)
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
%                                     It is relevant only in case parm(3)=1
%                    parm(3)      =1 ,  use Newton algorithm.
%                                    =2 , use Chebyshev algorithm.
%                    parm(4)      If parm(3)=2 then the interval which
%                                        includes the eigenvalues is [parm(4) parm(5)] 
%                    parm(5)      see parm(5).               
%
%   output
%                   u  =           the solution vectors
%
n=length(r);
tpoints=tpoints(:);
if parm(3)==1
     u=newton(r,tpoints,matvec,fun,m,tol,restarts,parm);
else
     u=chebalg(r,tpoints,matvec,fun,m,tol,parm);
end
