function y = wavelet(w, scale)
    % %Constants needed to compute the Morlet Wavelet
    w_0 = 6.0;
    c_sigma = (1 + exp(w_0.^2) - 2 * exp(-0.75*w_0.^2)).^(-0.5);
    sqrt_pi = pi.^-0.25;
    k_sigma = exp(-0.5 * w_0.^2);
    
    %Morlet Equation.
    y = c_sigma * sqrt_pi * ...
        ( exp(-0.5*(w_0 - scale * w).^2) - k_sigma * exp(-0.5* scale * w.^2));
end