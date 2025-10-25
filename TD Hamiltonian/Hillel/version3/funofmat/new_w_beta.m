function   [w,beta]=new_w_beta(w,beta,ro,ro0)
fac=ro/ro0;
p=1;
for j=length(w):-1:1
     w(j,:)=w(j,:)*p;
     p=p*fac;
     beta=beta/fac;
end