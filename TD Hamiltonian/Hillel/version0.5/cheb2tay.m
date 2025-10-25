function c=cheb2tay(a,t)
%the routine does the transformation  from Chebyshev expansion to Taylor expansion while 
% x is in the interval [0,t];
mu=2/t;
m=length(a);
alfa=zeros(m,1); beta=zeros(m,1);gama=zeros(m,1);
alfa(1)=1;alfa(2)=0;
beta(1)=-1;beta(2)=mu;
c=a(1)*alfa+a(2)*beta;
ju=[1:m]';
fac=2;
for k=1:m-2
    gama(1)=-2*beta(1)-alfa(1);
    gama(2:k+2)=2*mu*(ju(1:k+1).*beta(1:k+1))-2*beta(2:k+2)-alfa(2:k+2);
    alfa=beta;
    beta=gama;
    c=c+a(k+2)*gama;
end
 fac=1;
 for k=1:m
     c(k,:)=c(k,:)/fac;
     fac=fac*k;
 end

