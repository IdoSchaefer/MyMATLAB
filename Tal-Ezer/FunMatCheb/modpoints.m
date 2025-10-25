function x=modpoints(alpha,beta,n)
s1=g(1,alpha,beta);s2=g(-1,alpha,beta);
a=.5*(s1-s2);
b=.5*(s1+s2);
for j=1:n
    y(j,1)=cos((j-1)*pi/(n-1));
    x(j,1)=(1/a)*(g(y(j),alpha,beta)-b);
end
