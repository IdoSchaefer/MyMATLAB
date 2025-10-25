% 15'th harmonic target filter marking:
fill(w(1:1001)/0.06, exp(-(w(1:1001)-0.9).^2/(2*0.01^2)), [1, 0.8, 1], 'LineStyle', 'none')
fill(w(1:1001)/0.06, -exp(-(w(1:1001)-0.9).^2/(2*0.01^2)), [1, 0.8, 1], 'LineStyle', 'none')