%function test
n=60;
alpha=.8;beta=.1;
s1=g(1,alpha,beta);s2=g(-1,alpha,beta);
a=.5*(s1-s2);
b=.5*(s1+s2);
for j=1:n
    y(j,1)=cos((j-1)*pi/(n-1));
    x(j,1)=(1/a)*(g(y(j),alpha,beta)-b);
    % This is (dx/dy)^(-1):
    gdrv(j,1)=(a/sqrt(alpha*beta))*sqrt((1-alpha*y(j))*(1+beta*y(j)));
    % This is (d^x/dy^2)/(dx/dy):
    gttgt(j,1)=.5*(2*alpha*beta*y(j)+alpha-beta)/(1-alpha*y(j))/(1+beta*y(j));
    f(j,1)=fun(x(j));
    ff(j,1)=fun(y(j,1));
end
% one possibility for computing the second derivative
% w=chebdifft(f,2)-gttgt.*chebdifft(f,1);
% ap2=(gdrv.^2).*w;

%another possibility for computing the second derivative
w=gdrv.*chebdifft(f,1);
ap2=gdrv.*chebdifft(w,1);

%computing the second derivative by standard Chebyshev
apchb=chebdifft(ff,2);
for j=1:n
    ex(j,1)=dfun2(x(j));
    exchb(j,1)=dfun2(y(j));
end
er=norm(ex-ap2,inf)
erchb=norm(exchb-apchb,inf)
plot(x,ex,'r*')
hold on
plot(x,ap2,'k+')
