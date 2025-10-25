function fieldw_con = fieldw0b(fieldw_unc, T, filterE)
% The function computes a constrained field spectrum, fieldw_con, with 0 boundaries in the
% time domain, from the unconstrained field spectrum fieldw_unc.
    Nt = length(fieldw_unc) -  1;
    dw = pi/T;
    vfilterE = zeros(1, Nt + 1);
    for wi = 1:(Nt + 1)
        vfilterE(wi) = filterE((wi-1)*dw);
    end
    dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilterE(1), vfilterE(2:Nt), 0.5*vfilterE(Nt + 1)])*dw;
    coswT = ones(1, Nt + 1);
    coswT(2:2:(Nt + 1)) = -1;
    dctfilterET = sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1), vfilterE(2:Nt).*coswT(2:Nt), 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw;    
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    fieldt_unc0 = sqrt(2/pi)*sum([0.5*fieldw_unc(1), fieldw_unc(2:Nt), 0.5*fieldw_unc(Nt + 1)])*dw;    
    fieldt_uncT = sqrt(2/pi)*...
        sum([0.5*fieldw_unc(1)*coswT(1), fieldw_unc(2:Nt).*coswT(2:Nt), 0.5*fieldw_unc(Nt + 1)*coswT(Nt + 1)])*dw;    
    lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
    lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;
    fieldw_con = fieldw_unc - vfilterE.*(lambda0 + lambdaT*coswT);
end
