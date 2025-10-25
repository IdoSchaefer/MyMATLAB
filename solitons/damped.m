function x1dx2d=damped(t, x, g, k)
x1=x(1);
x2=x(2);
x1d=x2;
x2d=-g.*x2-k.*x1;
x1dx2d=[x1d;x2d];