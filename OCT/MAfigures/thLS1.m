% A 3LS example, with smaller penalty factor.
miu = [0 1 1;
       1 0 1;
       1 1 0];
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv, niter, mallniterc, J1, maxgrad, weight] = OCfMLSrlx([1;0;0], [1 1.9 3], [-1 5], miu,...
    @(w) exp(-20*(w-1).^2), @(w) 50*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), 0.1, 100, 0.05, 5, 5, 1e-3);
[fieldt, fieldw, psi, evmiut, evmiuw, relE, conv1, niter1, mallniterc, J1, maxgrad, weight] = OCfMLSrlx([1;0;0], [1 1.9 3], [-1 5], miu,...
    fieldw, @(w) 50*exp(-20*(w-1).^2), @(w) exp(-20*(w-2).^2), weight, 100, 0.05, 5, 5, 1e-3);
[evmiut0, evmiuw0, J10, psi0] = evmiuH0MLS(1/sqrt(2)*[1;0;1], diag([1 1.9 3]), miu, @(w) exp(-20*(w-2).^2), 100, 0.05);
[evmiuallnt, evmiuallnw, J1alln] = evmiuns(psi, miu, 100, @(w) exp(-20*(w-2).^2));
convtot=[conv, conv1(2:end)]; %size: 375
% J1=2.4292e+001,  J10=2.3989e+001
% J1>J10