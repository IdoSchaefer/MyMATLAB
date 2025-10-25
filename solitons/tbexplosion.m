clf
syms Fa Fb Fm Va Vb Vm r real;
getVa = input('Enter the potential function of the left chain in the terms of r: ', 's');
Va = eval(getVa);
Fa = -diff(Va);
getVb = input('Enter the potential function of the right chain in the terms of r: ', 's');
Vb = eval(getVb);
Fb = -diff(Vb);
getVm = input('Enter the potential function of the molecule''s bond in the terms of r: ', 's');
Vm = eval(getVm);
Fm = -diff(Vm);
N = input('Enter the number of the molecules: ');
%mdqsn = input('Do you want to insert mass defects?(y, n)?', 's');
%if mdqsn ~= 'y' && mdqsn ~= 'Y'
    ma = input('Enter the left bodie''s mass:');
    mb = input('Enter the right bodie''s mass:');
%else
%    sm = input('Enter the standard mass: ');
%    m = ones(1, N)*sm;
%    md = input('Enter the first kind of mass defect:');
%    nmd = input('Enter the number of the body with the mass defect (0 to continue): ');
%    while md>0 && nmd<=N && nmd>=1            
%        while nmd<=N && nmd>=1
%            m(nmd) = md;
%            nmd = input('Enter the number of the body with the mass defect (0 to continue): ');
%        end
%        md = input('Enter another kind of mass defect (0 to continue): '); 
%    end
%end
tfinal = input('Enter the time interval:');
solveFa = subs(solve(Fa));
r0 = min(solveFa);
r0mqsn = input('Choose: set the internal standard distance as the equilibrium distance of the molecule''s potential, \n or as your choice(e, c)?', 's');
if r0mqsn == 'c'
    r0m = input('Enter the standard distance:');
else
    solveFm = subs(solve(Fm));
    r0m = min(solveFm);
end
initxa = 0:r0:(N-1)*r0;
initxb = initxa + r0m;
initva = zeros(1, N);
initvb = zeros(1, N);
initvb(1) = input('Enter the initial v of the first right body: ');
sol = ode45(@fixedcexp2b, [0 tfinal], [initxa initxb initva initvb], [], Fa, Fb, Fm, r, ma, mb, N);
dt = 0.01;
t = 0:dt:tfinal;
x = deval(sol, t);
toptimei = fix(tfinal/dt + 1);
timei = 1:toptimei;
allravalues = x(2:N, timei) - x(1:(N-1), timei);
allrbvalues = x((N+2):(2*N), timei) - x((N+1):(2*N-1), timei);
allrvalues = allrbvalues;
maxra = max(max(allravalues));
maxrb = max(max(allrbvalues));
minra = min(min(allravalues));
minrb = min(min(allrbvalues));
tstep = 0.05;
tistep = 5;
hold on
i = 1:N;
plot(t, x(i,:))
plot(t, x(N+i,:))
xlabel('t (arbitrary units)')
ylabel('x (arbitrary units)')