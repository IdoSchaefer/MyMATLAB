function   [w,beta]=new_w_beta(w,beta,ro,ro0)
fac=ro/ro0;
p=1;
[s1,s2]=size(w);
for j=s1:-1:1
     w(j,:)=w(j,:)*p;
     p=p*fac;
     beta=beta/fac;
end