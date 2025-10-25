function result = torang4p(p1, p2, p3, p4)
% returns the torsion angle between the spaces specified by p1, p2, p3,
% and p2, p3, p4.
    v1 = p1 - p2;
    vax = p3 - p2;
    v2 = p4 - p3;
    if size(v1, 1) == 1
        v1 = v1';
        vax = vax';
        v2 = v2';
    end
% old version:
%    plane1 = [v1, vax];
%    plane2 = [vax, v2];    
%    result = subspace(plane1, plane2);
    result = angle2v(cross(cross(vax, v1), vax), cross(cross(vax, v2), vax));
end