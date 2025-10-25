function ww=leja(zz,w,k,ireal)
% The comments are mine (Ido).
% zz: The new candidate points.
% w: The old points
% k: The number of new points to be returned.
m=length(zz);
mu=length(w);
z=zz(1:m);
p=ones(m,1);
for j=1:mu
   p=(abs(z-w(j)).^(1/(mu+m))).*p;
end
j=mu+1;
while j <= mu+k
   [mx,jmx]=max(p);
   w(j,1)=z(jmx);
   if imag(z(jmx))==0 | ireal==0
       z(jmx:m-1)=z(jmx+1:m);
       p(jmx:m-1)=p(jmx+1:m);
       m=m-1;
       z=z(1:m,1); p=p(1:m,1);
       p=(abs(z-w(j)).^(1/(mu+m))).*p;
       j=j+1;
   else
       w(j+1,1)=conj(w(j,1));
       z(jmx:m-2)=z(jmx+2:m);
       p(jmx:m-2)=p(jmx+2:m);     
       m=m-2;
       z=z(1:m,1); p=p(1:m,1);
       p=(abs((z-w(j)).*(z-w(j+1))).^(1/(mu+m))).*p;
       j=j+2;
   end
end
 ww(:,1)=w(mu+1:mu+k);

   
