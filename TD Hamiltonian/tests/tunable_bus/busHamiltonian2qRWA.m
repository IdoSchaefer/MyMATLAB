% This is the same as busHamiltonian2qI, but in the framework of the RWA.
% The energy units are g1=2*pi*100MHz
omega1 = 5889.9/100;
omegac0 = 7445/100;
Dc0 = omegac0 - omega1;
D2 = 5031.1/100 - omega1;
alphac = 230/100;
alpha1 = 324/100;
alpha2 = 235/100;
g1 = 1;
g2 = 71.4/100;
n = (0:3).';
N_HO = spdiags(n, 0, 4, 4);
Manhar = spdiags([0; 0; -1; -3], 0, 4, 4);
a4 = spdiags(sqrt(n), 1, 4, 4);
adag4 = spdiags(sqrt(n + 1), -1, 4, 4);
%Mcoup = spdiags(sqrt([(n + 1), n]), [-1, 1], 4, 4);
I4 = speye(4);
H0I = multi_kron({Dc0*N_HO + alphac*Manhar, I4, I4}) + multi_kron({I4, alpha1*Manhar, I4})... 
    + multi_kron({I4, I4, D2*N_HO + alpha2*Manhar})...
    + g1*(multi_kron({adag4, a4, I4}) + multi_kron({a4, adag4, I4})) + g2*(multi_kron({adag4, I4, a4}) + multi_kron({a4, I4, adag4}));
% We shall utilize the symmetry of the Hamiltonian - it preserves the
% number of excitations. Since the columns of the unitary
% trasformation have a well defined number of excitations, we
% can use the relevant subspace only for each column. 
[H0excitations, Nsingle, Ndouble, singlei, doublei] = qubit_excitationsH(H0I, [4,4]);
T = 0.1142*100*2*pi;
theta = -0.108;
delta = 0.153;
sigma = 8.3e-3*100*2*pi;
omega_phi = 850.6/100;
envelope = @(t) exp(-(t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma))-T/2).^2/(2*sigma^2));
control = @(t) omegac0*(1 - sqrt(cos(pi*(theta + delta*envelope(t).*cos(omega_phi*t)))));
%control = @(t) omegac0*(1 - sqrt(cos(pi*envelope(t).*(theta + delta*cos(omega_phi*t)))));
Ncoupler = multi_kron({N_HO, I4, I4});
Ncoupler_excitations = qubit_excitationsH(Ncoupler, [4,4]);
GopRWA = @(u, t, v) -1i*(H0excitations*v - control(t)*(Ncoupler_excitations*v));
Gdiff_opRWA = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_excitations*u1));
Ntotal = 1 + 2*Nsingle + Ndouble;
u0RWA = zeros(Ntotal, 1);
indices2q = [1, 2, Nsingle + 3, 2*Nsingle + 3];
u0RWA(indices2q) = 1;
Usqrt_iswap_2q = sparse([1, 0,          0,          0;
                         0, 1/sqrt(2),  1i/sqrt(2), 0;
                         0, 1i/sqrt(2), 1/sqrt(2),  0;
                         0, 0,          0,          1]);
URWAsqrt_iswap = sparse(Ntotal, 1);
URWAsqrt_iswap([1, 2:3, Nsingle + (2:3), 2*Nsingle + 3]) = Usqrt_iswap_2q([1, 6:7, 10:11, 16]);
tgrid = 0:T/1e3:T;
options = SGdefault_op;
data = SGdata(options);
tic
[URWA, mniterRWA, matvecsRWA, est_errorsRWA, historyRWA] = SemiGlobal1(GopRWA, Gdiff_opRWA, 0, [], [-30*1i, 20*1i], u0RWA, tgrid, 2e3, 9, 2, 1e-5, options, data);
toc
% Restoring the full double-qubit columns of the unitary transformation:
restore_vec = [Nex_is_even; Nex_is_odd; Nex_is_odd; Nex_is_even];
% This takes local invariance into account:
fidelityRWA = sum(abs(URWA(:,end).*URWAsqrt_iswap(:)))^2/16;
infidelityRWA = 1 - fidelityRWA;