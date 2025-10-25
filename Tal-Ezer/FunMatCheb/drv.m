function  ap=drv(f)
global x g gtag
yd = chebdifft(f,1);
ap=g.*yd;
