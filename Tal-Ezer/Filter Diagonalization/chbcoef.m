function c=chbcoef(R,t)
m=floor(sqrt(3*R*t));
for j=1:m
    y(j)=cos((2*j-1)*pi/(2*m));
    z(j)=R*(y(j)-1)/2;
    f(j)=exp(R*t*z(j));
end
c=sqrt(2/m)*dct(f);
c(1)=c(1)/sqrt(2);
