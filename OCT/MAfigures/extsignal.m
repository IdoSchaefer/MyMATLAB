x=0:0.01:2*pi;
y=sin(x);
xext= -4*pi:0.01:6*pi;
yext= [sin(xext(1:floor(2*pi/0.01))), -sin(xext(ceil(2*pi/0.01):floor(4*pi/0.01))),...
    sin(xext(ceil(4*pi/0.01):floor(6*pi/0.01))), -sin(xext(ceil(6*pi/0.01):floor(8*pi/0.01))), sin(xext(ceil(8*pi/0.01):ceil(10*pi/0.01)))];
z=ones(1, 201)*2*pi;
figure
plot(x,y, 'k')
%here adjust width of the line
hold on
plot(xext, yext, 'g')
legend
plot(z, -1:0.01:1)
plot(-z, -1:0.01:1)