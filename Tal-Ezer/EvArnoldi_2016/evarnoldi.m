function [v,eg,er,flag]=evarnoldi(v,matvec,eta,n,k,m,tol,cycles,ireal)
% evarnoldi computes k eigenvalues of a matrix A defined by the routine matvec
%
% input:
%          v         - initial vector.
%          matvec - matvec is a routine written by the user. The head of the routine
%                        should be -     function u=matvec(v)
%                        where u,v satisfy  u=A*v
%          eta  -     if eta=-1  then the k smallest (in absolute value) eigenvalues are computed
%                       if eta=1  then the k largest (in absolute value) eigenvalues are computed
%                       if eta=alpha where alpha is not 1 or -1  then the k  eigenvalues, closest to alpha, are computed                           
%             n  -     the dimension of A.
%             k  -     the number of eigenvalues computed.
%            m   -    the dimension of the Krylov space.
%            tol -     acuuracy tolerance
%           cycles - maximal number of restarts.
%           ireal    -  if ireal=1, the matrix is real
%                    -   if ireal =2 , the matrix is complex
%
% output:
%           v - the computed eigenvectors
%           eg - the computed eigenvalues
%           er - the computed error
%           matvecs - number of matrix-vector multiplications performed.
%           flag -  if flag=0, the eigenvalues were computed to the desired accuracy.
%                   if flag=1, the eigenvalues were not computed to the desired accuracy.
%
zz=[];
flag=1;
v=v/norm(v);
h=[];
kp=floor(m/2);
for j=1:cycles
   [v,h]=lancz(v,h,matvec,n,m,j);
%    [v,h]=Earnoldi(v,h,matvec,m);
   [z,eg,q,p]=Eleja_points(h,kp,eta,ireal);
   [mp1,m]=size(h);
   er=h(m+1,m)*abs(q(m,1:k));
   res=max(er);
   if res < tol
        v=v(:,1:m)*q;
        flag=0;
        return
    end
    y=qpol(h,z); y(m+1)=0;
    if ireal==1
        y=real(y);
    end
    h(:,m+1)=zeros(m+1,1);
    [vs,h,ip]=sarnoldi(h,y);
    v=v*vs;
    if ip==1
        [q,dd]=eig(h);
        eg=diag(dd);
        v=v*q;
        flag=0;
        return
    end
 end
