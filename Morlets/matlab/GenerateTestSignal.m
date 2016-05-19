function signal = GenerateTestSignal(samples, FS, FREQ)

dt = 1.0./FS;
fsig = FREQ/FS;
dw = 2 * pi * fsig;
w0 = 0.01;
one_peri = 1./fsig;

signal = zeros(samples, 1);

signal(200:400) = sin( ((200:400) - 200)*dw + w0);
signal(1000:(1000+2*one_peri) ) = sin( ((1000:1000+2*one_peri) - 1000)*dw + w0);
signal(2000: (2000+3*one_peri) ) = sin( ((2000:2000+3*one_peri) - 2000)*dw + w0);
end