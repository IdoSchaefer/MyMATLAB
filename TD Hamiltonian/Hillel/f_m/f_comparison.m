function [rel_errors, allk] = f_comparison(m)
    rel_errors = zeros(1, 101);
    allk = zeros(1, 101);
    for zi = 1:101
        y_direct = my_direct_f(1i*(zi - 1)*0.1, m);
        [y_series, k] = f(1i*(zi - 1)*0.1, m);
        rel_errors(zi) = abs(y_direct - y_series)./abs(y_series);
        allk(zi) = k;
    end
    figure
    plot(0:0.1:10, log10(rel_errors))
    figure
    plot(0:0.1:10, allk)
end