% 11LS example 
miu = diag(ones(1, 10), 1) + diag(ones(1, 10), -1);
miu(1, 11)=1;
miu(11, 1)=1;
psi0=zeros(11, 1);
psi0(1)=1;
E0=[1 2.1 3 3.9 5 6.1 7 8.1 9 9.9 11];
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfMLSrlx(psi0, E0, [-1 13], miu,...
    @(w) rectanglefun(w, 0, 1.3), @(w) 50*rectanglefun(w, 0, 1.3), @(w) rectanglefun(w, 9.9, 10.1), 1, 100, 0.05, 9, 9, 1e-3);
