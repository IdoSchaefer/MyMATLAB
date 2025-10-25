function [u,tvec,iter,flag]=chb_tdrhs_exact1(u0,H,t,mcheb,tol,maxit)
global tvec
global mmm
mmm=mcheb;
n=length(u0);
flag=0;
iter=0;
iu=sqrt(-1);
tvec=tpoints(t,mcheb);
for j=1:mcheb
    u(:,j)=expm(H*tvec(j))*u0;
end
 %return
% for j=1:mcheb
%     u(:,j)=u0;
% end
ub=u;
fs=rhs1(u,tvec);
for  iter0=1:maxit
    a=dvd_vec(tvec,fs);
    s=convrt(tvec,a);
%     a= cheb(fs);
%     s=chbcf_convrt(a,t);
    v=u0;
    tj=ones(mcheb,1);
    for j=1:mcheb
        u(:,j)=u0;
    end
    for j=1:mcheb-1
        v=(H*v+s(:,j))/j;
        tj=tj.*tvec;
        u=u+v*tj';
    end
    v=(H*v+s(:,mcheb))/mcheb;
    for j=1:mcheb
         ww(:,j)=funmat(v,H,tvec(j),'fun');
    end
    u=u+ww;
    fs=rhs1(u,tvec);
     iter=iter+mcheb+1;
     er=norm(u-ub)/norm(ub);
     [iter0 norm(u) er];
     if er < tol
         return
     end
     ub=u;
     iter0;
end


    