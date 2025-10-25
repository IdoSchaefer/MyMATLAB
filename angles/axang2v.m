function result = axang2v(v1, v2, vax)
% calculates the angle between v1 and v2, arround the axis vector vax.
    vp1 = cross(vax, cross(v1, vax));
    vp2 = cross(vax, cross(v2, vax));
    result = angle2v(vp1, vp2);
end    