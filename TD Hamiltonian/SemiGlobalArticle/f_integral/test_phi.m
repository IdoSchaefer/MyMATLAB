function test_phi
n=30;
[x,wx]=gaussj(n,0,0);
z=[-100:-1]';
for m=1:20
    vm=(m/2)*((x+1)./2).^(m-1);
    for j=1:100
        f(:,j)=fexp(x,vm,z(j));
    end
    sol=f'*wx;
    hold on
    plot(z,sol)
end
