function GenerateTestSignal(samples, amplitude,streach, phase, noise)

signal = amplitude * cos( (1:samples) * streach + phase);

if noise > 0
    signal = awgn(signal,noise);
end

csvwrite('signal.csv', signal);
end