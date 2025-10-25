function z=real_leja(zset,k)
m=length(zset);
mmj=m;
zv=ones(mmj,1);
[a1,imax]=max(abs(zset));
j=0;
while j < k && mmj > 0
   j=j+1;
   z(j,1)=zset(imax);
   if imag(z(j)) ==0
       zset=zset([1:imax-1 imax+1:mmj]);
       zv=zv([1:imax-1 imax+1:mmj]);
       zv=abs(z(j)-zset).*zv;
       mmj=mmj-1;
   else
       zset=zset([1:imax-1 imax+2:mmj]);
       zv=zv([1:imax-1 imax+2:mmj]);
       z(j+1,1)=conj(z(j));
       zv=abs((z(j)-zset).*(z(j+1)-zset)).*zv;
       mmj=mmj-2;
       j=j+1;
   end
   if mmj==0
       return
   end
   [amax,imax]=max(zv);
end
