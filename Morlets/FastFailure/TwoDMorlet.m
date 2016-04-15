clc; close all; clear all;

%This should be a sample about 10 seconds long, sampled at 800 Hz
t = 10;
fs = 800;
samples = t * fs;

x = zeros(samples, 1);

%Add a 20 Hz cosine half way through the signal
cosine = cos(20*(1:20))';

x(2000:2019) = cosine;
x(4000:4019) = cosine;
x(6000:6019) = cosine;

%Generate the Morlet wavelet
% [psi, xMor] = morlet(-4, 4, 20);


%Convolute the Morlet wavelet with the signal
y = zeros(samples, 50);
scale = 1:50;
scale = 2.^scale;

for i = 1: length(scale)
    for j = 0: (samples/20 - 1)
        for k = 1:20
            y ( 20*j + k, i) = x(20*j + k) * Morlet(20*j + k, 800, scale(i)) ;
        end
    end
end