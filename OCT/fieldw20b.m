function fieldw_con = fieldw20b(fieldw, vfilterE, dw)
    Nt = length(fieldw) - 1;
    dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilterE(1); vfilterE(2:Nt); 0.5*vfilterE(Nt + 1)])*dw;
    coswT = ones(Nt + 1, 1);
    coswT(2:2:(Nt + 1)) = -1;
    dctfilterET = sqrt(2/pi)*sum([0.5*vfilterE(1)*coswT(1); vfilterE(2:Nt).*coswT(2:Nt); 0.5*vfilterE(Nt + 1)*coswT(Nt + 1)])*dw;
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    fieldt_unc0 = sqrt(2/pi)*sum([0.5*fieldw(1); fieldw(2:Nt); 0.5*fieldw(Nt + 1)])*dw;    
    fieldt_uncT = sqrt(2/pi)*...
            sum([0.5*fieldw(1)*coswT(1); fieldw(2:Nt).*coswT(2:Nt); 0.5*fieldw(Nt + 1)*coswT(Nt + 1)])*dw;    
    lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
    lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;
    fieldw_con = fieldw - vfilterE.*(lambda0 + lambdaT*coswT);
end