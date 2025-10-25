function newv = rotv(v, ang, axisv)
    szv = size(v);
    if szv(1) == 1
        v = v';
    end
    rotM = [cos(ang) -sin(ang)  0;
            sin(ang)  cos(ang)  0;
            0         0         1];
    if nargin > 2
        szaxv = size(axisv);
        if szaxv(1) == 1
           axisv = axisv';
        end
        if axisv(1) == 0 && axisv(2) == 0
            if axisv(3) < 0
                rotM(1, 2) = -rotM(1, 2);
                rotM(2, 1) = -rotM(2, 1);
            end
        else
            moveax = zeros(3, 3);
            moveax(:, 3) = axisv/norm(axisv);
            moveax(1:2, 2) = [-moveax(2, 3); moveax(1, 3)];
            moveax(:, 2) = moveax(:, 2)/norm(moveax(:, 2));
            moveax(:, 1) = cross(moveax(:, 2), moveax(:, 3));
            rotM = Mtrans(rotM, moveax);
        end
    end
    newv = rotM*v;
    if szv(1) == 1
        newv = newv';
    end
end
        