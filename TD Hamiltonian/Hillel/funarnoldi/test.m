function test
global A
global matvecs
matvecs=0;
n=128;
A=100*matrix(n);
% load fHH HH
% [n,n]=size(HH);
% A=HH;
% A=mmread('orani678.mtx');
% [n,n]=size(A);
b=rand(n,1);
save fb b
load fb b
t=1; jp=1:n;maxit=1000;tol=1.e-5;m=10;
[x,flag]=fun_of_matrix(b,'matvec','fun',t,m,tol,maxit,jp);
for j=1:length(t)
    ex(:,j)=expm(A*t(j))*b;
    er(j,1)=norm(ex(:,j)-x(:,j));
end
er
matvecs
