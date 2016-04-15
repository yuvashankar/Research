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
[psi, xMor] = morlet(-4, 4, 20);
% plot(xMor, psi)

%Convolute the Morlet wavelet with the signal
y = zeros(samples, 50);

plot(y(1950:2020))
