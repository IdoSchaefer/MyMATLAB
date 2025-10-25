function eqx = Vequilibrium(Ffun, field)
    Nt = length(field);
    eqx = zeros(1, Nt);
    defaultop = optimset('fsolve');
    options = optimset(defaultop, 'Display', 'off');
%     r0 = fsolve(eqfun, 0, options);
%     r1 = fsolve(@(x) eqfun(x) - 1, r0, options);
    eqx(1) = fsolve(@(x) Ffun(x) + field(1), 0, options);
    for ti = 2:Nt
        eqx(ti) = fsolve(@(x) Ffun(x) + field(ti), eqx(ti - 1), options);
    end
end