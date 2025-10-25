function u=rk4(u,t,dt)
global x
u0=u;
iu=sqrt(-1);
fs=ff(t)*(x.*u);
u=(iu/2)*(fourdifft(u,2)-(x.*x).*u);
u=u-iu*fs;
k1=dt*u;

u1=u0+k1/2;
fs=ff(t+dt/2)*(x.*u1);
u=(iu/2)*(fourdifft(u1,2)-(x.*x).*u1);
u=u-iu*fs;
k2=dt*u;


u1=u0+k2/2;
fs=ff(t+dt/2)*(x.*u1);
u=(iu/2)*(fourdifft(u1,2)-(x.*x).*u1);
u=u-iu*fs;
k3=dt*u;

u1=u0+k3;
fs=ff(t+dt)*(x.*u1);
u=(iu/2)*(fourdifft(u1,2)-(x.*x).*u1);
u=u-iu*fs;
k4=dt*u;


u=u0+(1/6)*(k1+2*(k2+k3)+k4);

