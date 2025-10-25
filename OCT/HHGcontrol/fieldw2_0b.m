function fieldw0b = fieldw2_0b(fieldw, vfilter)
% The function computes a field spectrum, fieldw0b, constrained to be with 0 boundaries in the
% time domain, from the unconstrained field spectrum fieldw. The new field
% lies in the spectral range specified by a filter function, represented by
% the vector vfilter
    Nw = length(fieldw) - 1;
    coswT = ones(1, Nw + 1);
    coswT(2:2:(Nw + 1)) = -1;
    fieldt_unc0 = sqrt(2/pi)*sum([0.5*fieldw(1), fieldw(2:Nw), 0.5*fieldw(Nw + 1)]);    
    fieldt_uncT = sqrt(2/pi)*...
    sum([0.5*fieldw(1)*coswT(1), fieldw(2:Nw).*coswT(2:Nw), 0.5*fieldw(Nw + 1)*coswT(Nw + 1)]);    
    dctfilterE0 = sqrt(2/pi)*sum([0.5*vfilter(1), vfilter(2:Nw), 0.5*vfilter(Nw + 1)]);
    dctfilterET = sqrt(2/pi)*sum([0.5*vfilter(1)*coswT(1), vfilter(2:Nw).*coswT(2:Nw), 0.5*vfilter(Nw + 1)*coswT(Nw + 1)]);    
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
    lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;
    fieldw0b = fieldw - vfilter.*(lambda0 + lambdaT*coswT);
end
