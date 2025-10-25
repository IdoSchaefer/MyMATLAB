function solx = frog(getDx , initx, t0tf, dt, params)
    sizeix = size(initx);
    if sizeix(1)==1
        initx = initx';
    end
    N = max(sizeix);
    initDx = getDx(initx, params);
    toptimei = fix((t0tf(2)-t0tf(1))/dt) + 1;
    solx = zeros(N, toptimei);
    solx(:, 1) = initx;
    solx(:, 2) = initx + initDx*dt;
    for ti = 2:(toptimei - 1)
        Dxt = getDx(solx(:, ti), params);
        solx(:, ti + 1) = solx(:, ti - 1) + Dxt*2*dt;
    end