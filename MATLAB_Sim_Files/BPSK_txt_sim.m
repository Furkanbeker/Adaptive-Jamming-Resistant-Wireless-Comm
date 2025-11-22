clc; clear all; close all;

% File Path Settings
inputPath = 'C:\Users\mfb36\Desktop\input.txt';
outputPath = 'C:\Users\mfb36\Desktop\output.txt';

% Read text from input.txt
fid = fopen(inputPath, 'r');
if fid == -1
    error('Unable to open input.txt. Check the file path.');
end

% ASCII to Bit Sequence Conversion
binData = dec2bin(fread(fid, '*char')', 8);
txBits = reshape(binData.' - '0', 1, []);
fclose(fid);

% BPSK Modulation
symbols = 2 * txBits - 1;

% Timing settings
br = 1e6; f = br; T = 1/br;
fs = 100 * f;
t = 0:1/fs:T-1/fs;
samplesPerSymbol = length(t);

% Carrier modulation
modulatedSignal = [];
for i = 1:length(txBits)
    s = symbols(i) * cos(2*pi*f*t);
    modulatedSignal = [modulatedSignal s];
end

% Signal Power Calculation
signal_power = mean(modulatedSignal.^2);

% Set Noise Power Manually 
noise_power = 0.1; 

% Noise Addition
noise = sqrt(noise_power) * randn(1, length(modulatedSignal));
noisySignal = modulatedSignal + noise;

% SNR Calculation (from given noise_power)
calculated_SNR_linear = signal_power / noise_power;
calculated_SNR_dB = 10 * log10(calculated_SNR_linear);

% Demodulation
rxBits = zeros(1, length(txBits));
for i = 1:length(txBits)
    idx_start = (i-1)*samplesPerSymbol + 1;
    idx_end = i*samplesPerSymbol;
    r = noisySignal(idx_start:idx_end);
    r_demod = r .* cos(2*pi*f*t);
    r_avg = mean(r_demod);
    rxBits(i) = r_avg > 0;
end

% BER Calculation
excessBits = mod(length(rxBits), 8);
if excessBits ~= 0
    rxBits = rxBits(1:end - excessBits);
    txBits = txBits(1:end - excessBits);
end
ber = sum(txBits ~= rxBits) / length(rxBits);

fprintf('Signal Power: %.4f\n', signal_power);
fprintf('Noise Power: %.4f\n', noise_power);
fprintf('Calculated SNR: %.2f dB\n', calculated_SNR_dB);
fprintf('Bit Error Rate (BER): %.4f\n', ber);


% Write received characters to output.txt
rxBitsMatrix = reshape(rxBits, 8, []).';
rxChars = char(bin2dec(num2str(rxBitsMatrix)))';
fid_out = fopen(outputPath, 'w');
fwrite(fid_out, rxChars);
fclose(fid_out);

% Plotting
signalTime = linspace(0, T*length(txBits), length(modulatedSignal));

figure(1);
subplot(4,1,1);
plot(signalTime, modulatedSignal, 'LineWidth', 1.2);
title('BPSK In-phase Component (I)');
ylabel('Amplitude'); grid on;

subplot(4,1,2);
plot(signalTime, zeros(size(modulatedSignal)), 'LineWidth', 1.2);
title('BPSK Quadrature Component (Q)');
ylabel('Amplitude'); grid on;

subplot(4,1,3);
plot(signalTime, noise, 'LineWidth', 1.2);
title('Channel Noise');
ylabel('Amplitude'); grid on;

subplot(4,1,4);
plot(signalTime, noisySignal, 'r', 'LineWidth', 1.2);
title('Noisy BPSK Signal');
xlabel('Time'); ylabel('Amplitude'); grid on;

figure(2);
subplot(2,1,1);
stem(txBits(1:min(100, length(txBits))), 'filled');
title('Transmitted Bits'); ylim([-0.5 1.5]); grid on;

subplot(2,1,2);
stem(rxBits(1:min(100, length(rxBits))), 'filled');
title('Received Bits (Decisions)'); ylim([-0.5 1.5]); grid on;

figure(3);
scatter(noisySignal(1:samplesPerSymbol:end), zeros(1,length(txBits)), 200, 'rx');
title('BPSK Constellation Diagram');
xlabel('In-phase'); ylabel('Quadrature');
xline(0, 'k'); yline(0, 'k');
xlim([-2 2]); ylim([-0.5 0.5]); axis square; grid on;

