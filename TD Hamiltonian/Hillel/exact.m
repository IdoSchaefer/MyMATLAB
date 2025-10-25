function u=exact(t)
global mu
global E0
global T
iu=sqrt(-1);
u(1,1)=cos(1/4*mu*E0*(t-T/2/pi*sin(2*pi*t/T)));
u(2,1)=-iu*sin(1/4*mu*E0*(t-T/2/pi*sin(2*pi*t/T)));