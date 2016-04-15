function y = Morlet (x, w0, scale)

    w02 = w0.^2;
    normal = 1./sqrt(scale);
    k = exp(-0.5*w02);
    rootPi = pi.^(-1/4);
    
    x = x*scale;
    
    y = normal * rootPi * exp(-0.5 * x * x) * cos( w0 * x - k);
end