function [allNt, allmv, aller, all_est_ers] = error_bus2qSG(Nt_ts, Nfm, minNt, Nsamp, Niter, Arnoldi)
    if nargin<5
        Niter = 1;
    end
    if nargin<6
        Arnoldi = false;
    end
    if ~Arnoldi
        ev_domain = [-528*1i, 1i];
    else
        ev_domain = [];
    end
% The energy units are g1=2*pi*100MHz
omegac0 = 7445/100;
omega1 = 5889.9/100;
omega2 = 5031.1/100;
alphac = 230/100;
alpha1 = 324/100;
alpha2 = 235/100;
g1 = 1;
g2 = 71.4/100;
n = (0:3).';
N_HO = spdiags(n, 0, 4, 4);
Manhar = spdiags([0; 0; -1; -3], 0, 4, 4);
Mcoup = spdiags(sqrt([(n + 1), n]), [-1, 1], 4, 4);
I4 = speye(4);
H0 = multi_kron({omegac0*N_HO + alphac*Manhar, I4, I4}) + multi_kron({I4, omega1*N_HO + alpha1*Manhar, I4})... 
    + multi_kron({I4, I4, omega2*N_HO + alpha2*Manhar})...
    + g1*multi_kron({Mcoup, Mcoup, I4}) + g2*multi_kron({Mcoup, I4, Mcoup});
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
H0even = H0(Nex_is_even, Nex_is_even);
H0odd = H0(Nex_is_odd, Nex_is_odd);
% The propagation Hamiltonian for the unitary transformation:
H0u = sparse(blkdiag(H0even, H0odd, H0odd, H0even));
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
Gop = @(u, t, v) -1i*(H0u*v - control(t)*(Ncoupler_u*v));
Gdiff_op = @(u1, t1, u2, t2) 1i*(control(t1) - control(t2)).*(Ncoupler_u*u1);
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
[allNt, allmv, aller, all_est_ers] = error_decaySG2(Gop_odd, Gdiff_op_even, 0, [], ev_domain, U2q0par(:, 2), Upar_ex2qT(33:64), [0, T], Nt_ts, Nfm, minNt, Nsamp, Niter);
end