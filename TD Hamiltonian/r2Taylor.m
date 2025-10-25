function Q = r2Taylor(x)
    NC = length(x);
    Q = zeros(NC);
    % first row:
    Q(1, 1) = 1;
    % all the rest of the rows:
    for ri = 2:NC
        Q(ri, 1) = -x(ri - 1)*Q(ri - 1, 1);
        for Taylori = 2:ri-1
            Q(ri, Taylori) = (Taylori - 1)*Q(ri - 1, Taylori - 1) - x(ri - 1)*Q(ri - 1, Taylori);
        end
        Q(ri, ri) = (ri - 1)*Q(ri - 1, ri - 1);
    end
end    