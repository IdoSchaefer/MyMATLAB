function options = options_OCf_qn(tolx, maxNiter)
% The function generates the options structure for the quasiNewton procedure
% for optimal control of harmonic generation.
%%%% An old version; I haven't used it for long. It has to be updated. Use
%%%% optionsOCqn.m instead.
    options.invHess0 = [];
    options.externalHessian = false;
    options.finvHess = @HessianBFGS;
    options.f_termination = @(dif_f, dif_x, fmin, gradmin, xmin) norm(dif_x)/norm(xmin)<=tolx;
    options.maxNiter = maxNiter;
    options.maxNiter_ls_1st = 20;
    options.maxNiter_ls = 10;
    options.ro = 0;
    options.sigma = 0.9;
    options.tau1 = 9;
    options.tau2 = 0.1;
    options.tau3 = 0.5;
    options.minimal_f = -Inf;
    options.Deltaf0 = [];
    options.plot = @(x, y) plot(x, -y);
end