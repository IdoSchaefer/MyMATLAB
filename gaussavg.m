function yavg = gaussavg(y, sigma, nsigma, dx)
    N = length(y);
    Lgaussian = nsigma*sigma;
    gaussian = exp(-0.5*((0:dx:Lgaussian)/sigma).^2);
    n0g = length(gaussian);
    gaussian = [gaussian(end:-1:2), gaussian];
    Ng = length(gaussian);
    for yi = 1:N
    end
end
        
        sumg = 