function ww=real_lejaf(zz,w,k)
[m,m2]=size(zz);
if m2 > 1
    m=m2;
    zz=conj(zz');
end
mu=length(w);
z=zz(1:m);
p=ones(m,1);
for j=1:mu
   p=(abs(z-w(j)).^(1.d0/(mu+m))).*p;
end
j=mu+1;
while j <= mu+k
   [mx,jmx]=max(p);
   w(j,1)=z(jmx);
   z(jmx:m-1)=z(jmx+1:m);
   p(jmx:m-1)=p(jmx+1:m);
   m=m-1;
   z=z(1:m,1); p=p(1:m,1);
   if imag(z(jmx))==0
       p=(abs(z-w(j)).^(1.d0/(mu+k))).*p;
       j=j+1;
   else
       w(j+1,1)=conj(w(j,1));
       p=(abs((z-w(j)).*(z-w(j+1))).^(1.d0/(mu+k))).*p;
       j=j+2;
   end
end
 ww(:,1)=w(mu+1:mu+k);

   
