clc; close all; clear all;
% % Initializing the constants
fs = 1000; %sampling rate
dt = 1./fs;
N = 3 * fs;
frequency = 1;
disc_freq = frequency/fs;
dw = 2 * pi * disc_freq;
phase_shift = 0.01;

% %create a 3 second cosine signal to test
test_input = 0:N;
t = test_input * dt;
test_input = cos(dw * test_input + phase_shift);

%Prepare the Daughter Wavelet array to fit the way Matlab's FFT works
deltaW =(2 * pi)/N;
% deltaW = 6000/N;
w = -N/2:N/2;
w =  deltaW * w;

% %Constants needed to compute the Morlet Wavelet
w_0 =  6.0;
scale = (2 * pi) * 1;
c_sigma = (1 + exp(-w_0.^2) - 2 * exp(-0.75*w_0.^2)).^(-0.5);
sqrt_pi = pi.^-0.25;
k_sigma = exp(-0.5 * w_0.^2);

%Morlet Equation.
y = c_sigma * sqrt_pi * ...
    ( exp(-0.5*(w_0 - scale * w).^2) - k_sigma * exp(-0.5 * scale * w.^2));


hold on
plot(w,y);
xlim([-pi pi]);
ylabel('omega (rad/s)');