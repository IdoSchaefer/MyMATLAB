function [z, ro]=leja_points1(h,zint,ro)
[mp1,m]=size(h);
zp=[];
for j=1:m
    eg=eig(h(1:j,1:j));
    zp=[zp;eg];
end
if ro == 0
    ro=rofun(zp);
end    
 z=leja_complex(zp,zint,m,ro);
%  z=leja_complex(zp,m);