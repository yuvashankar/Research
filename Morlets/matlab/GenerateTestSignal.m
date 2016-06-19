function GenerateTestSignal(samples, amplitude,streach, phase, noise)

signal = amplitude * cos( (1:samples) * streach + phase);

if noise > 0
    signal = awgn(signal,noise);
end
signal = signal';
% csvwrite('signal.csv', signal);
save ( 'signal.txt', 'signal', '-ascii');
end