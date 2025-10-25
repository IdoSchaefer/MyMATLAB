function c=cheb_coef(a,b,fun,n,m,t)
for i=1:n
    tet=cos((2*i-1)*pi/2/n);
    x(i)=.5*((b-a)*tet+a+b);
    f(i,1)=feval(fun,x(i),m,t);
end
c=sqrt(2/n)*dct(f);
c(1)=c(1)/sqrt(2);
j=1;
iu=sqrt(-1);
for i=1:n
    c(i)=j*c(i);
    j=j*iu;
end
c=real(c);