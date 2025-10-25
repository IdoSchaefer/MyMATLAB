% This is the same as busHamiltonian2q, but in the interaction picture of
% H_1 = omega1*(N_c + N_1 + N_2). As expected, the convergence is
% slower, but the function of matrix approximation requires less terms.
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
%H0u = kron(I4, H0);
% We shall utilize the symmetry of the Hamiltonian - it preserves the
% parity of the number of excitations. Since the columns of the unitary
% trasformation have a well defined parity of the number of excitations, we
% can use the relevant subspace only for each column. 
% The number of excitations in each entry of the direct product space:
Nex = Nexcitations_kron([4, 4], (1:64).');
% Classification into even and odd number of excitations:
Nex_is_even = (mod(Nex, 2) == 0);
Nex_is_odd = ~Nex_is_even;
% Extracting the even and odd number of excitation terms from the full
% Hamiltonian:
H0Ieven = H0I(Nex_is_even, Nex_is_even);
H0Iodd = H0I(Nex_is_odd, Nex_is_odd);
% The propagation Hamiltonian for the unitary transformation:
H0Iu = sparse(blkdiag(H0Ieven, H0Iodd, H0Iodd, H0Ieven));
T = 0.1142*100*2*pi;
theta = -0.108;
delta = 0.153;
sigma = 8.3e-3*100*2*pi;
omega_phi = 850.6/100;
envelope = @(t) exp(-(t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma))-T/2).^2/(2*sigma^2));
control = @(t) omegac0*(1 - sqrt(cos(pi*(theta + delta*envelope(t).*cos(omega_phi*t)))));
%control = @(t) omegac0*(1 - sqrt(cos(pi*envelope(t).*(theta + delta*cos(omega_phi*t)))));
Ncoupler = multi_kron({N_HO, I4, I4});
Ncoupler_even = Ncoupler(Nex_is_even, Nex_is_even);
% We utilize the fact that the even part of Ncoupler is identical to the
% odd one, and thus the odd part needn't be computed.
Ncoupler_u = kron(I4, Ncoupler_even);
%Ncoupler_u = kron(I4, Ncoupler);
b_a1 = multi_kron({a4, a4, I4});
b_a2 = multi_kron({a4, I4, a4});
bdag_a1dag = multi_kron({adag4, adag4, I4});
bdag_a2dag = multi_kron({adag4, I4, adag4});
b_a1even = b_a1(Nex_is_even, Nex_is_even);
b_a1odd = b_a1(Nex_is_odd, Nex_is_odd);
b_a2even = b_a2(Nex_is_even, Nex_is_even);
b_a2odd = b_a2(Nex_is_odd, Nex_is_odd);
bdag_a1dag_even = bdag_a1dag(Nex_is_even, Nex_is_even);
bdag_a1dag_odd = bdag_a1dag(Nex_is_odd, Nex_is_odd);
bdag_a2dag_even = bdag_a2dag(Nex_is_even, Nex_is_even);
bdag_a2dag_odd = bdag_a2dag(Nex_is_odd, Nex_is_odd);
b_a1u = sparse(blkdiag(b_a1even, b_a1odd, b_a1odd, b_a1even));
b_a2u = sparse(blkdiag(b_a2even, b_a2odd, b_a2odd, b_a2even));
bdag_a1dag_u = sparse(blkdiag(bdag_a1dag_even, bdag_a1dag_odd, bdag_a1dag_odd, bdag_a1dag_even));
bdag_a2dag_u = sparse(blkdiag(bdag_a2dag_even, bdag_a2dag_odd, bdag_a2dag_odd, bdag_a2dag_even));
GopI = @(u, t, v) -1i*(H0Iu*v - control(t)*(Ncoupler_u*v) +...
    exp(1i*2*omega1*t)*(g1*(bdag_a1dag_u*v) + g2*(bdag_a2dag_u*v)) + exp(-1i*2*omega1*t)*(g1*(b_a1u*v) + g2*(b_a2u*v)));
Gdiff_opI = @(u1, t1, u2, t2) -1i*((control(t2) - control(t1)).*(Ncoupler_u*u1) +...
    (exp(1i*2*omega1*t1) - exp(1i*2*omega1*t2)).*(g1*bdag_a1dag_u*u1 + g2*bdag_a2dag_u*u1) +...
    (exp(-1i*2*omega1*t1) - exp(-1i*2*omega1*t2)).*(g1*b_a1u*u1 + g2*b_a2u*u1));
U2q0 = generateU2q0full([4, 4, 4]);
% The active entries of U2q0:
U2q0par = [U2q0(Nex_is_even, 1),U2q0(Nex_is_odd, 2), U2q0(Nex_is_odd, 3), U2q0(Nex_is_even, 4)];
Usqrt_iswap_2q = sparse([1, 0,          0,          0;
                         0, 1/sqrt(2),  1i/sqrt(2), 0;
                         0, 1i/sqrt(2), 1/sqrt(2),  0;
                         0, 0,          0,          1]);
Usqrt_iswap = U2q0*Usqrt_iswap_2q;
tgrid = 0:T/1e3:T;
options = SGdefault_op;
data = SGdata(options);
tic
[UIpar, mniterI, matvecsI, est_errorsI, historyI] = SemiGlobal1(GopI, Gdiff_opI, 0, [], [-50*1i, 50*1i], U2q0par(:), tgrid, 10e3, 9, 2, 1e-5, options, data);
toc
% Restoring the full double-qubit columns of the unitary transformation:
restore_vec = [Nex_is_even; Nex_is_odd; Nex_is_odd; Nex_is_even];
UI = zeros(256, 1001);
UI(restore_vec, :) = UIpar;
% This takes local invariance into account:
fidelityI = sum(abs(UI(:,end).*Usqrt_iswap(:)))^2/16;
infidelityI = 1 - fidelityI;