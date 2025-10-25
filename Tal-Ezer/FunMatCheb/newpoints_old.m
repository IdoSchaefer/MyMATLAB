function znew=newpoints(z,h,ireal,method)
[mp1,m]=size(h);
zz=[];
if method==1
    for j=2:m
        e=eig(h(1:j,1:j));
        zz=[zz;e];
    end
elseif method==2
    for j=2:m-1
        hh=h(1:j,1:j)*h(1:j,1:j);
        hh(:,j)=hh(:,j)+h(j+1,j)*h(1:j,j+1);
        e=-sqrt(eig(hh));
        zz=[zz;e];
    end
end
% zz=cfilter(zz);
if ireal==1
    kz=length(zz);
    k=0;
    j=0;
    while j < kz
        j=j+1;
        k=k+1;
        zzz(k,1)=zz(j);
        if imag(zz(j)) ~=0
            j=j+1;
        end
    end
    zz=zzz;
end
znew=real_lejaf(zz,z,m);
if imag(sum(znew))~=0 & ireal==1
    znew(m+1)=conj(znew(m));
  % znew=znew(1:m-1);
end
    
