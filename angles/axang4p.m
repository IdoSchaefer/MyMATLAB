function result = axang4p(axp1, axp2, p1, p2)
% calculates the angle between p1-axp1, p2-axp1, arround the axis vector axp2-axp1.
    v1 = p1 - axp1;
    v2 = p2 - axp1;
    vax = axp2 - axp1;
    vp1 = cross(vax, cross(v1, vax));
    vp2 = cross(vax, cross(v2, vax));
    result = angle2v(vp1, vp2);
end    