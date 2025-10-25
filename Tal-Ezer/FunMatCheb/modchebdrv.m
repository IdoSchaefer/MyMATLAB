function w=modchebdrv(u,alpha,beta)
n=length(u);
s1=g(1,alpha,beta);s2=g(-1,alpha,beta);
a=.5*(s1-s2);
b=.5*(s1+s2);
for j=1:n
    y(j,1)=cos((j-1)*pi/(n-1));
    x(j,1)=(1/a)*(g(y(j),alpha,beta)-b);
    gdrv(j,1)=(a/sqrt(alpha*beta))*sqrt((1-alpha*y(j))*(1+beta*y(j)));
    gttgt(j,1)=.5*(2*alpha*beta*y(j)+alpha-beta)/(1-alpha*y(j))/(1+beta*y(j));
end
% one possibility for computing the second derivative
w=chebdifft(u,2)-gttgt.*chebdifft(u,1);
w=(gdrv.^2).*w;

%another possibility for computing the second derivative
% w=gdrv.*chebdifft(u,1);
% w=gdrv.*chebdifft(w,1);
