function y=fun(z,t)
global mglobal
m=mglobal;
nz=length(z);
nt=length(t);
for iz=1:nz
    for it=1:nt
        zt=z(iz)*t(it);
        if abs(zt) < m+1
            y(iz,it)=1;
            r=1;
            k=m;
            while abs(r)/abs(y(iz,it))> 1.e-13
                k=k+1;
                r=r*zt/k;
                y(iz,it)=y(iz,it)+r;
            end
        else
            y(iz,it)=exp(zt);
            for j=1:m
                y(iz,it)=j*(y(iz,it)-1)/zt;
            end
        end
        y(iz,it)=y(iz,it)*t(it)^m;
    end
end

