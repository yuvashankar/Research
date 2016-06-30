function y = FourierMorlet(w, w0, scale, n)

    heavisideMatrix = heaviside(w);

    scaleMatrix(1:n) = scale;
%     scaleMatrix = scaleMatrix';
    w0Matrix(1:n) = w0;
%     w0Matrix = w0Matrix';


    exponent = (scaleMatrix * w - w0Matrix) * (scaleMatrix * w - w0Matrix) * heavisideMatrix;
    norm = sqrt(scaleMatrix * w(1)) * (pi.^-0.25) * sqrt(n);

    y = norm * exp(exponent);
    y = y * heavisideMatrix;
end
