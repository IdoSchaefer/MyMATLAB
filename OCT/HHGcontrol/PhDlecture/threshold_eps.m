Nfield = 100;
threshold_field = zeros(1, Nfield);
xthreshold_field = zeros(1, Nfield);
for fieldi = 1:Nfield
    field = fieldi*1e-3;
    xthreshold_field(fieldi) = fzero(@(x) x/(1 + x^2)^(3/2) - field, [0.5, 1e2]);
    threshold_field(fieldi) = 1 - 1/sqrt(1 + xthreshold_field(fieldi)^2) - xthreshold_field(fieldi)*field;
end