niter1ab = niter1a+niter1b
conv1ab = [conv1a conv1b(2:end)];
figure
plot(0:niter1ab, conv1ab)
hold on
plot([7 7], [0 7e-3])
xlabel('number of iterations', 'interpreter', 'latex')
ylabel('$J$','interpreter', 'latex')