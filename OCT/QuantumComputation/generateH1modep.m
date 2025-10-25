load pcontrol2qubits
a7s = spdiags(sqrt((0:6).'), 1, 7, 7);
adag7s = spdiags(sqrt((1:7).'), -1, 7, 7);
Manhar = spdiags(Enlg, 0, 7, 7);
I5s = speye(5);
H05 = multi_kron({I5s, Manhar(1:5, 1:5), I5s}) + multi_kron({I5s, I5s, Manhar(1:5,1:5)});
p5 = 1i*(adagq(1:5,1:5) - aq(1:5,1:5));
Hd1 = multi_kron({I5s, p5, I5s});
Hd2 = multi_kron({I5s, I5s, p5});
adag7s = spdiags(sqrt((1:7).'), -1, 7, 7);
Hcplus = multi_kron({a7s(1:5, 1:5), adagq(1:5, 1:5), I5s}) + multi_kron({a7s(1:5, 1:5), I5s, adagq(1:5, 1:5)});
Hcminus = Hcplus';
I4s = speye(4);
H05u = kron(I4s, H05);
Hd1u = kron(I4s,Hd1);
Hd2u = kron(I4s,Hd2);
Hcplusu = kron(I4s,Hcplus);
Hcminusu = kron(I4s,Hcminus);
Hoperations.psi = @(psi, params, v) H05u*v - params(1)*(Hd1u*v) - params(2)*(Hd2u*v) + params(3)*(Hcplusu*v) + conj(params(3))*(Hcminusu*v);
Hoperations.chi = Hoperations.psi;
Hoperations.diff_psi = @(psi1, params1, psi2, params2) (params2(1) - params1(1, :)).*(Hd1u*psi1) + (params2(2) - params1(2, :)).*(Hd2u*psi1)...
    + (params1(3, :) - params2(3)).*(Hcplusu*psi1) + conj(params1(3, :) - params2(3)).*(Hcminusu*psi1);
Hoperations.diff_chi = Hoperations.diff_psi;
fcouplingOp = {@(v) Hd1u*v; @(v) Hd2u*v};
gs5 = sparse(5 ,1);
gs5(1) = 1;
fundamental5 = sparse(5, 1);
fundamental5(2) = 1;
U0 = [multi_kron({gs5, gs5, gs5}), multi_kron({gs5, gs5, fundamental5}), multi_kron({gs5, fundamental5, gs5})...
    multi_kron({gs5, fundamental5, fundamental5})];
Utarget = [multi_kron({gs5, gs5, gs5}), 1i*multi_kron({gs5, fundamental5, gs5}), 1i*multi_kron({gs5, gs5, fundamental5}), ...
    multi_kron({gs5, fundamental5, fundamental5})];
ones5 = ones(5,1);
phi4 = sparse(5, 1);
phi4(5) = 1;
penalforbv = double(multi_kron({ones5, phi4, ones5}) | multi_kron({ones5, ones5, phi4}));
penalforbvu = kron(ones(4,1),penalforbv);
fJforb1 = get_fJforb_penal(penalforbvu);
options = optionsOCqn(1e-4, 1e4);
filterE = @(w) 0.5*(1-tanh(2*(w-13))).*heaviside(15-w);
options.f_max_alpha = get_f_max_alphaOCf_multiE(5, 13*pi/2e3, 13*pi, filterE, 2);