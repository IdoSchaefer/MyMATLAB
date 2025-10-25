function tvec=tpoints(t,m)
for j=1:m
 %   s=cos((2*j-1)*pi/(2*m));
    s=-cos((j-1)*pi/(m-1));
    tvec(j,1)=t*(1+s)/2;
end
