function v=chebfun(v,A,thalf,m,a,b)
% the routine approximates w=f(A,t)v by w_m = P_m(A,t)v
% where P_m(z,t) is Chebyshev expansion polynomial which 
% approximates f(z,t) in [a,b] which is the inteval which 
% includes all the eigenvalues of A.
% input:
%          v = input vectors
%          A = input operator
%          t = vector of time points
%          fun = input function
%          m  = maximal degree of Chebyshev polynomial
%          [a,b]= inpit interval of eigenvalues
%
% output: 
%           w=output vectors
global thalfglobal
thalfglobal=thalf;
v(:,1)=v;
v(:,2)=(1/(b-a))*(2*feval(A,v(:,1))-(a+b)*v(:,1));
for j=3:m
    v(:,j)=(2/(b-a))*(2*feval(A,v(:,j-1))-(a+b)*v(:,j-1))-v(:,j-2);
end

