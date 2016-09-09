
FS = 1000.0;
dt = 1.0/FS;
FREQ = 19.0;
n = 3000;

signal = GenerateTestSignal(n, FS, FREQ);

fftSignal = fft(signal);

k = [-floor(n/2 - 1):1:floor(n/2)]';
kplus = kplus * (2*pi/(n + dt));

y = FourierMorlet(k, 5.0, 22.0, n);
