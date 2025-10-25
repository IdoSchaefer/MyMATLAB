function znew=newpoints(z,h,ireal)
% [mp1,m]=size(h);
% zz=[];
% for j=2:m
%     e=eig(h(1:j,1:j));
%     zz=[zz;e];
% end
% znew=leja(zz,z,m,ireal);
% if imag(sum(znew))~=0 & ireal==1
%    znew(m+1)=conj(znew(m));
% %    znew=znew(1:m-1);
% end
% return
[mp1,m]=size(h);
zz=[];
for j=2:m-1
    hh=h(1:j,1:j)*h(1:j,1:j);
    hh(:,j)=hh(:,j)+h(j+1,j)*h(1:j,j+1);
    e0=sqrt(eig(hh));
    e1=-abs(real(e0)); e2=-abs(imag(e0));
    ee=complex(e1,e2);
    k=0;
    for i=1:j
        k=k+1;
        e(k,1)=ee(i);
        if real(ee(i))==0 
            k=k+1;
            e(k,1)=conj(ee(i));
        end
    end
    zz=[zz;e];
end
znew=leja(zz,z,m,ireal);
if imag(sum(znew))~=0 & ireal==1
    znew(m+1)=conj(znew(m));
  %   znew=znew(1:m-1);
end
    
