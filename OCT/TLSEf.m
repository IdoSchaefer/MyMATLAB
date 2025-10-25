function Dpsi = TLSEf(t, psi, E, field, T)
% Linear interpolation between time steps:
%     midti = t/T*(length(field)-1) + 1;
%     floorti = floor(midti);
%     ceilti = ceil(midti);
%     weight = midti-floorti;
%     field_t = field(floorti)*(1-weight) + field(ceilti)*weight;
    dt = T/(length(field) - 1);
    field_t = RK4interp(field, 0:dt:T, t);
    Ht = diag(E) + [0         -conj(field_t);
                    -field_t  0             ];
    Dpsi = -1i*Ht*psi;
end
    