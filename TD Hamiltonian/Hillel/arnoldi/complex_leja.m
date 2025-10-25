function z=complex_leja(zset,k)
m=length(zset);
mmj=m;
zv=ones(mmj,1);
[a1,imax]=max(abs(zset));
j=0;
while j < k && mmj > 0
   j=j+1;
   z(j,1)=zset(imax);
   zset=zset([1:imax-1 imax+1:mmj]);
   zv=zv([1:imax-1 imax+1:mmj]);
  zv=abs(z(j)-zset).*zv;
  mmj=mmj-1;
  if mmj==0
       return
   end
   [amax,imax]=max(zv);
end
