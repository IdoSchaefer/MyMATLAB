%function test4
%n=30;
n=129;
alpha=.99;
 beta=.99;
s1=g(1,alpha,beta);s2=g(-1,alpha,beta);
a=.5*(s1-s2);
b=.5*(s1+s2);
for j=1:n
    y(j,1)=cos((j-1)*pi/(n-1));
    gdrv(j,1)=(a/sqrt(alpha*beta))*sqrt((1-alpha*y(j))*(1+beta*y(j)));
    gttgt(j,1)=.5*(2*alpha*beta*y(j)+alpha-beta)/(1-alpha*y(j))/(1+beta*y(j));
end
u=zeros(n,1);
% My comment:
% Computation of the second derivative matrix A. Each column is computed by
% the opperation of the second derivative operator on a delta function,
% represented by u. The 1'st column is computed by the operation on the 0
% vector, probably to enforce boundary conditions.
for j=1:n
    u(j)=1;
    u(1)=0;
%     w=gdrv.*chebdifft(u,1);
%     A(:,j)=gdrv.*chebdifft(w,1);
%     A(n,j)=0;
% My comment - computation of the operation of the second derivative operator on u.  
    w=chebdifft(u,2)-gttgt.*chebdifft(u,1);
    A(:,j)=(gdrv.^2).*w;
    A(n,j)=0;
    u(j)=0;
end
eg=sort(eig(A))
plot(real(eg),imag(eg),'r*')
axis equal
