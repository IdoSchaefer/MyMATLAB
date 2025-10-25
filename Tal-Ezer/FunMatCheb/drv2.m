function  ap=drv2(f)
global x g gtag
yd = chebdifft(f,1);
yd2=chebdifft(f,2);
ap=(g.*g).*yd2+yd.*gtag;
