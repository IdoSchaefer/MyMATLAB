% The energy units are g1=2*pi*100MHz
omegac0 = 7445/100;
omega1 = 5889.9/100;
omega2 = 5031.1/100;
omega3 = 4350/100;
alphac = 230/100;
alpha1 = 324/100;
alpha2 = 235/100;
alpha3 = 1;
g1 = 1;
g2 = 71.4/100;
g3 = 85/100;
n = (0:3).';
N_HO = spdiags(n, 0, 4, 4);
Manhar = spdiags([0; 0; -1; -3], 0, 4, 4);
Mcoup = spdiags(sqrt([(n + 1), n]), [-1, 1], 4, 4);
I4 = speye(4);
% The order of the sub-systems in the direct product:
% Tunable bus coupler, transmon 3, transmon 1, transmon 2
H03q = multi_kron({omegac0*N_HO + alphac*Manhar, I4, I4, I4}) + multi_kron({I4, I4, omega1*N_HO + alpha1*Manhar, I4})... 
    + multi_kron({I4, I4, I4, omega2*N_HO + alpha2*Manhar}) + multi_kron({I4, omega3*N_HO + alpha3*Manhar, I4, I4})...
    + g1*multi_kron({Mcoup, I4, Mcoup, I4}) + g2*multi_kron({Mcoup, I4, I4, Mcoup}) + g3*multi_kron({Mcoup, Mcoup, I4, I4});
% We shall utilize the symmetry of the Hamiltonian - it preserves the
% parity of the number of excitations. Since the columns of the unitary
% trasformation have a well defined parity of the number of excitations, we
% can use the relevant subspace only for each column. 
% The number of excitations in each entry of the direct product space:
Nex3q = Nexcitations_kron([4, 4, 4], (1:256).');
% Classification into even and odd number of excitations:
Nex_is_even3q = (mod(Nex3q, 2) == 0);
Nex_is_odd3q = ~Nex_is_even3q;
% Extracting the even and odd number of excitation terms from the full
% Hamiltonian:
H03q_even = H03q(Nex_is_even3q, Nex_is_even3q);
H03q_odd = H03q(Nex_is_odd3q, Nex_is_odd3q);
% The propagation Hamiltonian for the unitary transformation:
H03q_u = sparse(blkdiag(H03q_even, H03q_odd, H03q_odd, H03q_even));
T = 0.1142*100*2*pi;
theta = -0.108;
delta = 0.153;
sigma = 8.3e-3*100*2*pi;
omega_phi = 850.6/100;
envelope = @(t) exp(-(t + 0.5*log(cosh(t-T+3*sigma)./cosh(t-3*sigma))-T/2).^2/(2*sigma^2));
control = @(t) omegac0*(1 - sqrt(cos(pi*(theta + delta*envelope(t).*cos(omega_phi*t)))));
%control = @(t) omegac0*(1 - sqrt(cos(pi*envelope(t).*(theta + delta*cos(omega_phi*t)))));
Ncoupler3q = multi_kron({N_HO, I4, I4, I4});
Ncoupler3q_even = Ncoupler3q(Nex_is_even3q, Nex_is_even3q);
% We utilize the fact that the even part of Ncoupler3q is identical to the
% odd one, and thus the odd part needn't be computed.
Ncoupler3q_u = kron(I4, Ncoupler3q_even);
Gop3q = @(u, t, v) -1i*(H03q_u*v - control(t)*(Ncoupler3q_u*v));
Gdiff_op3q = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler3q_u*u1);
U2q0_3q = generateU2q0full([4, 4, 4, 4]);
% The active entries of U2q0:
U2q0par3q = [U2q0_3q(Nex_is_even3q, 1) ,U2q0_3q(Nex_is_odd3q, 2), U2q0_3q(Nex_is_odd3q, 3), U2q0_3q(Nex_is_even3q, 4)];
Usqrt_iswap_2q = sparse([1, 0,          0,          0;
                         0, 1/sqrt(2),  1i/sqrt(2), 0;
                         0, 1i/sqrt(2), 1/sqrt(2),  0;
                         0, 0,          0,          1]);
Usqrt_iswap = U2q0_3q*Usqrt_iswap_2q;
tgrid = 0:T/1e3:T;
options = SGdefault_op;
data = SGdata(options);
tic
[Upar3q, mniter3q, matvecs3q, est_errors3q, history3q] = SemiGlobal1(Gop3q, Gdiff_op3q, 0, [], [-655*1i, 1i], U2q0par3q(:), tgrid, 1e4, 9, 5, 1e-5, options, data);
toc
% Restoring the full double-qubit columns of the unitary transformation:
restore_vec3q = [Nex_is_even3q; Nex_is_odd3q; Nex_is_odd3q; Nex_is_even3q];
U3q = zeros(1024, 1001);
U3q(restore_vec3q, :) = Upar3q;
% This takes local invariance into account:
fidelity3q = sum(abs(U3q(:,end).*Usqrt_iswap(:)))^2/16;
infidelity3q = 1 - fidelity3q;