fieldwg2 = 3*(exp(-(w-w(19)).^2/(2*0.01^2)).*sin((w-w(19))*pi/0.015) + exp(-(w-w(10)).^2/(2*0.005^2)).*sin((w-w(10))*pi/0.0075));
fieldtg2 = dctI(fieldwg2)/dctfactor;