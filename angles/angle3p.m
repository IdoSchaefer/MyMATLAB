function ang = angle3p(p1, p2, p3)
    v1 = p1 - p2;
    v2 = p3 - p2;
    szv1 = size(v1);
    if szv1(1) > 1
        v1 = v1';
    end
    szv2 = size(v2);
    if szv2(1) == 1
        v2 = v2';
    end    
    ang = acos(v1/norm(v1)*v2/norm(v2));
end