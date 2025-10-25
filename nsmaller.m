function vi = nsmaller(v, x)
    vi = 1;
    while v(vi) > x
        vi = vi + 1;
    end
end
