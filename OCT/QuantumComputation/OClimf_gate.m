function [fieldt, fieldw, psi, relE, conv, niter, mallniterc, Jterms, maxgrad, alpha, invHess] = OClimf_gate(psi0, target, limits,...
    Hoperations, Hdiff_matvecs, fcouplingOp, Edomain, fguess, filterE, options, T, dt, Nt_ts, Nfm, tolprop)
% The function solves an OCT procedure for the optimization of the
% final-time expectation value problem, with limitations on the frequency band of
% the driving field.
% A quasi-Newton optimization procedure is applied.
%%% Input: %%%
% psi0: The initial state vector
% Vf: The stationary potential; a functon handle of the form @(x), where
% the input is the x grid. Alternatively, Vf can be a column vector of the
% potential values in x.
% m: The mass (in HHG, the electron mass, which is 1 in atomic units)
% xdomain: The boundaries of the x domain, where xdomain = [min_x, max_x]
% mu_x: The dipole function mu evaluated at the x grid points; a column 
% vector of the dimension of the x grid; mu_x is used in the Hamiltonian for the
% interaction term with the external field, -field(t)*mu_x. In HHG, insert 
% here the numerical x vector which its derivative decays to 0 near the 
% absorbing boundaries.
% a0: the stationary acceleration, evaluated at the x grid points; a
% column vector of the dimension of the x grid 
% penalnormf: The penalty function on the ionization; a function handle of
% the form @(sq_norm_psi), where the input is <psi(T)|psi(T)>
% Dpenalnormf: The derivative of the function represented by penalnormf; a
% function handle of the form @(sq_norm_psi), where the input is <psi(T)|psi(T)>
% fguess: The guess field function in the frequency domain; a function
% handle of the form @(w), where w is a column vector which represents the
% omega grid. Aternatively, fguess can be a column vector of the guess
% function evaluated at the omega grid points. If the guess does not
% satisfy the 0 boundary conditions of the temporal profile, the program
% constructs a constrained field with the 0 boundary conditions from the
% unconstrained field.
% filterE: The scaled filter function of the driving electric field spectrum; a
% function handle of the form @(w), where w is a column vector which
% represents the omega grid
% filtera: The filter function of the emission spectrum, represented by the
% stationary acceleration expectation spectrum; a function handle of the
% form @(w), where w is a column vector which represents the omega grid
% options: The options structure of the BFGS search procedure (see
% quasiNewton.m and default_op_qn.m); if is empty (substituted with []), the default is set by
% the procedure optionsOCqn.m, using the tol and maxNiter input parameters.
% T: The final time
% dt: The time-step
% Nt_ts: The number of Chebyshev points in the internal grid in each
% time-step, using the semi-global propagator.
% Nkr: The dimension of the Krylov approximation in the semi-global propagator
% tolprop: The tolerance parameter of the semi-global propagator
%%% Output: %%%
% fieldt: The optimized time-dependent electric field, evaluated at the
% points of the main time-grid, t = 0:dt:T (not at the intermediate 
% Chebyshev point structure of the semi-global propagation); a row vector
% fieldw: The optimized electric field in the omega domain; a row vector
% psi: The state in all time-points of the process under the influence of the optimized
% field; evaluated at the main time-grid, where the different time-points
% are represented by separate columns.
% evat: The time-dependent expectation value of the stationary
% acceleration, evaluated at the main time-grid; a row vector
% evaw: evat in the omega domain; a row vector
% evmut: The time-dependent expectation value of the dipole, defined by
% mu_x; a row vector
% evmuw: evmut in the omega domain; a row vector
% relE: The relative difference of the field from the previous iteration at the
% end of the optimization process
% conv: The convergence history, i.e. J as a function of the iteration number;
% a row vector of dimension niter+1, where the 1'st entry represents the
% guess field
% niter: The total number of iterations
% mallniterc: The mean number of semi-global iterations in the whole
% optimization process; should be close to 1 for efficiency.
% Jterms: A structure which contains the values of the different components
% of J. It has the fields: Jmax, Jenergy, Jion.
% maxgrad: The infinity norm of the gradient at the end of the optimization
% process
% alpha: The alpha parameter at the end of the optimization process (see
% quasiNewton.m)
% invHess: The inverse Hessian at the end of the optimization process

    % The number of time-points:
    Nt = round(T/dt);
    % The distance between adjacent points in the omega grid:
    dw = pi/T;
    % Integration weights of the boundary including omega grid:
    integw = [dw/2; ones(Nt - 1, 1)*dw; dw/2];
    % Additional factor for conversions between the time and frequency grid
    % by the DCT transformation:
    dctfactor = T/(sqrt(Nt*pi));
    dim = size(psi0);
    % The dimension of the Hilbert space:
    N_Hilbert = dim(1);
    % The index of the last point in the semi-global time-grid, where the
    % internal Chebyshev points are also counted:
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    % The index of the last point in the semi-global time-grid, including
    % the test-point for error computation:
    allt_tp_lasti = Nt*Nt_ts + 1;
    % The time indices of the propagation grid excluding the test-points:
    indices_excluding_tp = 1:allt_tp_lasti;
    indices_excluding_tp(Nt_ts:Nt_ts:Nt*Nt_ts) = [];
%     is_not_testp = true(1, allt_tp_lasti);
%     % The entries of the propagation grid, without the test points, are
%     % true:
%     is_not_testp(Nt_ts:Nt_ts:Nt*Nt_ts) = false;
    % The indices rearranged for the backward propagation:
    backwardi = zeros(1, allt_tp_lasti);
    backwardi(indices_excluding_tp) = fliplr(indices_excluding_tp);
    backwardi(Nt_ts:Nt_ts:Nt*Nt_ts) = Nt*Nt_ts:-Nt_ts:Nt_ts;
    % The sum of the mean number of semi-global iterations of all 
    % propagations performed during the optimization process (required for
    % the computation of mallniterc):
    summniterc = 0;
    % The number of control terms:
    Ncontrols = length(fcouplingOp);
    % The time-dependent field in the semi-global propagation grid (including the internal
    % Chebyshev points in each time-step):
    allfield = zeros(Ncontrols, allt_tp_lasti);
    % The different terms in the functional J:
    Jmax = 0;
    Jenergy = 0;
%     % The term -imag[<chi(t)|mu|psi(t)>] evaluated at the semi-global propagation
%     % grid:
%     chi_coup_psi = zeros(allt_lasti, Ncontrols);
    % The DCT of -imag[<chi(t)|mu|psi(t)>]:
    dct_chi_coup_psi = zeros(Nt + 1, Ncontrols);
    % The reversed Chebyshev points in each time-step in the [-1, 1] domain:
    tcheb = -cos(((1:Nt_ts).' - 1)*pi/(Nt_ts-1));
    % The internal Chebyshev points of each time-step in the [0, dt] domain:
    t_ts = 0.5*(tcheb+1)*dt;
    tmidi = round(Nt_ts/2);
    % The test time-point for time-expansion error computation (see the SemiGlobal programs):
    test_tpoint = (t_ts(tmidi) + t_ts(tmidi + 1))/2;
    % The time-step internal structure including the test-point:
    t_ts_tp = [t_ts(1:(Nt_ts - 1)); test_tpoint];
    % The test-point for the backward propagation (coincides with the
    % test-point of the forward propagation, but is defined from the
    % opposite propagation direction):
    tp_backward = -(t_ts(tmidi) + t_ts(tmidi - 1))/2;
    % The omega grid:
    w = (0:pi/T:pi/dt).';
    % The vector of the driving field filter-function, evaluated at w:
    vfilterE = filterE(w);
    % The indices of the w which participate in the optimization, where
    % the driving field spectrum does not attain negligible values (nz is
    % non-zero):
    iEnz = find(vfilterE>=eps*max(vfilterE));
    NEnz = length(iEnz);
%    [iEnz, ~, vfilterEnz] = find(vfilterE);
%    iEnz = find(vfilterE>=tolprop*max(vfilterE)/10);
    % vfilterE in the points which participate in the optimization:
    vfilterEnz = vfilterE(iEnz);
    % The integration weights of the internal grid Chebyshev points in the
    % semi-global prapagation time grid:
    igweights = chebweights(Nt_ts, 1).';
    igweights = [2*igweights(1); igweights(2:(Nt_ts - 1))];
    % The integration weights of the omega grid points which participate in
    % the optimization:
    integwnz = integw(iEnz);
    % The inverse cosine transform of the scaled filter function, evaluated
    % at t=0:
    dctfilterE0 = sqrt(2/pi)*sum(vfilterEnz.*integwnz);
    % Constructing a vector of cos(w*T) evaluated at all w values:
    coswT = ones(Nt + 1, 1);
    coswT(2:2:(Nt + 1)) = -1;
    % The inverse cosine transform of the scaled filter function, evaluated
    % at t=T:
    dctfilterET = sqrt(2/pi)*sum(vfilterEnz.*coswT(iEnz).*integwnz);
    % The determinant of the M matrix in the paper:
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    % The driving field spectrum in all w:
    fieldw = zeros(Nt + 1, Ncontrols);
    if isa(fguess, 'function_handle')
    % If fguess is a function handle:
        fieldw(iEnz, :) = repmat(fguess(w(iEnz)), 1, Ncontrols);
    elseif isa(fguess, 'cell')
    % If fguess is a cell of function handles:
        for ifguess = 1:Ncontrols
            fieldw(iEnz, ifguess) = fguess{ifguess}(w(iEnz));
        end
    else
    % If fguess is a vector or a matrix:
        fieldw(iEnz, :) = fguess(:, iEnz).';
    end
    % The driving field spectrum evaluated at the omega values which participate in the
    % optimization:
    fieldwnz = fieldw(iEnz, :);
    % If the guess field doesn't satisfy the temporal boundary conditions, a new
    % guess field is constructed which satisfies the boundary conditions:
    if sqrt(2/pi)*abs(sum(fieldwnz.*integwnz)) > tolprop*1e-2 | sqrt(2/pi)*abs(sum(fieldwnz.*coswT(iEnz).*integwnz)) > tolprop*1e-2
        fieldwnz = getfieldwcon(fieldwnz);
    end
    % psi at the semi-global propagation grid:
    allpsi = zeros(N_Hilbert, allt_tp_lasti);
    % nprop counts the number of propagations:
    nprop = 0;
    if isempty(options)
        % Setting the default options:
        options = optionsOCqn(1e-4, 1e3);
    end
    if isempty(options.invHess0)
        % Initializing the approximated Hessian to the default: 
        options.invHess0 = diag(repmat(vfilterEnz./(2*integwnz), Ncontrols, 1));
    end
    % fieldwnz is vectorized for the use of the quasiNewton procedure, which
    % is defined for vectors only:
    fieldwnz = fieldwnz(:);
    if options.maxNiter > 0
        % Performing an optimization process of fieldwnz:
        [fieldwnz, ~, minusgrad, niter, ~, ~, dif_fieldw, minus_conv, alpha, invHess] = quasiNewton(@Jeval, fieldwnz, options);
    else
        % If options.maxNiter == 0, the program returns the output variables
        % for the input field:
        [minus_conv, minusgrad] = Jeval(fieldwnz);
        niter = 0;
        dif_fieldw = 0;
        alpha = [];
        invHess = options.invHess0;
    end
    % Preparing the output variables:
    fieldw = fieldw.';
    maxgrad = max(abs(minusgrad));
    conv = -minus_conv;
    relE = norm(dif_fieldw)/norm(fieldwnz);
    niter
    nprop
    psi = allpsi(:, 1:Nt_ts:allt_tp_lasti);
    fieldt = allfield(:, 1:Nt_ts:allt_tp_lasti);
    mallniterc = summniterc/nprop;
    Jterms.Jmax = Jmax;
    Jterms.Jenergy = Jenergy;
    beep
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad] = Jeval(fieldwnz)
    % The optimization function
    % Input:
    % fieldwnz: The driving field spectrum at the omega points which
    % participate in the optimization
    % Output:
    % minusJ: The objective of the optimization, which is -J
    % minusgrad: The gradient of the optimization, which is minus the
    % gradient of J w.r.t. fieldwnz
        fieldwnzM = reshape(fieldwnz, NEnz, Ncontrols);
        fieldw(iEnz, :) = fieldwnzM;
        % Obtaining the the temporal driving field from its spectrum,
        % evaluated at the semi-global propagation grid:
        for icontrol = 1:Ncontrols
            allfield(icontrol, :) = dctIintgrid1(fieldw(:, icontrol).', T, t_ts_tp.')/dctfactor;            
        end
        % The forward propagation of psi:
        [~, mniterc, ~, max_errors_psi, history_psi] = SemiGlobalHparams(Hoperations.psi, Hoperations.diff_psi, Hdiff_matvecs,...
            allfield, [], Edomain, psi0, [0, T], Nt, Nt_ts, Nfm, tolprop, 10, 16, test_tpoint, false);
        allpsi(:, indices_excluding_tp) = history_psi.U;
        allpsi(:, Nt_ts:Nt_ts:Nt*Nt_ts) = history_psi.Utestp;
        summniterc = summniterc + mniterc;
        psiT = allpsi(:, allt_tp_lasti);
        % The final condition for the chi backward propagation:
        [chiT, Jmax] = chis_gate(psiT, target, limits);
        % The backward propagation of chi:
        [~, mniterc, ~, max_errors_chi, history_chi] = SemiGlobalHparams(Hoperations.chi, Hoperations.diff_chi, Hdiff_matvecs,...
            allfield(:, backwardi), [], Edomain, chiT, [T, 0], Nt, Nt_ts, Nfm, tolprop, 10, 16, tp_backward, false);
        allchi = history_chi.U;
        % Note that allchi doesn't include the test-points, unlike allpsi.
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        % Computing the objective -J:
        Jenergy_fun = -(fieldwnzM.^2./vfilterEnz).*integwnz;
        Jenergy = sum(sum(Jenergy_fun));
        minusJ = -(Jmax + Jenergy);
        % Computing the -gradient.
        if Ncontrols == 1
            % For a single control term:
            %  Obtaining -imag[<chi(t)|mu|psi(t)>] at the propagation grid:
            chi_coup_psi = -imag(sum(conj(allchi(:, allt_lasti:-1:1)).*fcouplingOp(allpsi(:, indices_excluding_tp)))).';
            % Computing chi_coup_psi transformed to the frequency domain,
            % utilizing the data in the entire propagation grid:
            dct_chi_coup_psi = dctfactor*dctIfrom_ig_sym1(chi_coup_psi, T, t_ts(1:(Nt_ts - 1)), igweights);
        else
            % For multiple control terms:
            for icontrol = 1:Ncontrols
                chi_coup_psi = -imag(sum(conj(allchi(:, allt_lasti:-1:1)).*fcouplingOp{icontrol}(allpsi(:, indices_excluding_tp)))).';
                dct_chi_coup_psi(:, icontrol) = dctfactor*dctIfrom_ig_sym1(chi_coup_psi, T, t_ts(1:(Nt_ts - 1)), igweights);
            end
        end
        % The unconstrained fieldw; a matrix of dimension NEnz*Ncontrols:
        newfieldw_unc_nzM = dct_chi_coup_psi(iEnz, :).*vfilterEnz;
        % Computing the Lagrange-multiplies of the temporal boundary
        % constraints from the unconstrained field. Note that if there are
        % multiple controls, then lambda0 and lambdaT become row vectros of the
        % dimension of Ncontrols:
        [lambda0, lambdaT] = fieldw_unc2lambda(newfieldw_unc_nzM);
        % The gradient arranged in a matrix of dimension NEnz*Ncontrols:
        minusgradM = 2*integwnz.*(fieldwnzM./vfilterEnz - dct_chi_coup_psi(iEnz, :) + lambda0 + coswT(iEnz)*lambdaT);
        % minusgradM is vectorized for the use of the quasiNewton
        % procedure:
        minusgrad = minusgradM(:);
    end

    function [lambda0, lambdaT] = fieldw_unc2lambda(fieldw_unc_nzM)
    % The function computes the Lagrange-multiplies of the boundary
    % constraints, given an unconstrained field in the frequency domain,
    % evaluated at the omega values which participate in the optimization.
        % Computing the unconstrained field at t=0, by the DCT of the
        % driving field specrum evaluated at t=0:
        fieldt_unc0 = sqrt(2/pi)*sum(fieldw_unc_nzM.*integwnz);
        % Computing the unconstrained field at t=T, by the DCT of the
        % driving field specrum evaluated at t=0:
        fieldt_uncT = sqrt(2/pi)*sum(fieldw_unc_nzM.*coswT(iEnz).*integwnz);
        lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
        lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;        
    end

    function fieldw_con = getfieldwcon(fieldw_unc_nz)
    % The function computes a constrained field spectrum, fieldw_con, with 0 boundaries in the
    % time domain, from the unconstrained field spectrum fieldw_unc.
        [lambda0, lambdaT] = fieldw_unc2lambda(fieldw_unc_nz);
        fieldw_con = fieldw_unc_nz - vfilterEnz.*(lambda0 + coswT(iEnz)*lambdaT);
    end

end