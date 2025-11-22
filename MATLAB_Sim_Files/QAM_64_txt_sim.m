clc; clear all; close all;

% File Path Settings
inputPath = 'C:\Users\mfb36\Desktop\input.txt';
outputPath = 'C:\Users\mfb36\Desktop\output.txt';

% Read text from input.txt
fid = fopen(inputPath, 'r');
if fid == -1
    error('input.txt file could not be opened.');
end
txt = fread(fid, '*char')';
fclose(fid);

% ASCII → Bit Conversion
binData = dec2bin(txt, 8);
txBits = reshape(binData.' - '0', 1, []);
numBits = length(txBits);
if mod(numBits,6) ~= 0
    txBits(end+1 : end+mod(6 - mod(numBits,6), 6)) = 0;
    numBits = length(txBits);
end

% 64QAM Modulation
I_levels = [-7 -5 -3 -1 1 3 5 7];
Q_levels = [-7 -5 -3 -1 1 3 5 7];
I_bin = de2bi(0:7, 3, 'left-msb');
Q_bin = de2bi(0:7, 3, 'left-msb');
symbols = zeros(1, numBits/6);

for i = 1:6:numBits
    b = txBits(i:i+5);
    I_idx = find(ismember(I_bin, b(1:3), 'rows'));
    Q_idx = find(ismember(Q_bin, b(4:6), 'rows'));
    I_val = I_levels(I_idx);
    Q_val = Q_levels(Q_idx);
    symbols((i+5)/6) = (I_val + 1j * Q_val) / sqrt(42);
end

% Signal and Noise Power Setup
signal_power = mean(abs(symbols).^2);
noise_power = 0.1;
snr_linear = signal_power / noise_power;
SNR_dB = 10 * log10(snr_linear);
noise_sym = sqrt(noise_power/2) * (randn(1,length(symbols)) + 1j*randn(1,length(symbols)));
rx_symbols = symbols + noise_sym;

% Demodulation
rxBits = zeros(1, numBits);
for i = 1:length(rx_symbols)
    I_val = real(rx_symbols(i)) * sqrt(42);
    Q_val = imag(rx_symbols(i)) * sqrt(42);
    [~, I_idx] = min(abs(I_val - I_levels));
    [~, Q_idx] = min(abs(Q_val - Q_levels));
    rxBits(6*i-5:6*i-3) = I_bin(I_idx,:);
    rxBits(6*i-2:6*i)   = Q_bin(Q_idx,:);
end

% Remove extra bits and calculate BER
excess = mod(length(rxBits), 8);
if excess ~= 0
    rxBits = rxBits(1:end-excess);
    txBits = txBits(1:end-excess);
end
ber = sum(rxBits ~= txBits) / length(rxBits);

% Display results
fprintf('Signal Power: %.4f\n', signal_power);
fprintf('Noise Power: %.4f\n', noise_power);
fprintf('SNR (dB): %.2f\n', SNR_dB);
fprintf('Bit Error Rate (BER): %.4f\n', ber);

% Bit → ASCII and Write to output.txt
rxBitsMatrix = reshape(rxBits, 8, []).';
rxChars = char(bin2dec(num2str(rxBitsMatrix)))';
fid_out = fopen(outputPath, 'w');
fwrite(fid_out, rxChars);
fclose(fid_out);

% Physical Signal Generation
br = 1e6; f = br; T = 1/br;
fs = 100 * f;
t = 0:1/fs:T-1/fs;
y_mod = []; y_in = []; y_qd = [];

for i = 1:length(symbols)
    I = real(symbols(i)) * sqrt(42);
    Q = imag(symbols(i)) * sqrt(42);
    in = I * cos(2*pi*f*t);
    qd = Q * sin(2*pi*f*t);
    y_in = [y_in in];
    y_qd = [y_qd qd];
    y_mod = [y_mod in + qd];
end

noise_time = sqrt(noise_power) * randn(1, length(y_mod));
y_noisy = y_mod + noise_time;
tt = linspace(0, T*length(symbols), length(y_mod));

% Plots
figure(1);
subplot(4,1,1);
plot(tt, y_in, 'LineWidth', 1.2);
title('64QAM In-phase Component (I)');
ylabel('Amplitude'); grid on;

subplot(4,1,2);
plot(tt, y_qd, 'LineWidth', 1.2);
title('64QAM Quadrature Component (Q)');
ylabel('Amplitude'); grid on;

subplot(4,1,3);
plot(tt, noise_time, 'LineWidth', 1.2);
title('Channel Noise');
ylabel('Amplitude'); grid on;

subplot(4,1,4);
plot(tt, y_noisy, 'r', 'LineWidth', 1.2);
title('Noisy 64QAM Signal');
xlabel('Time'); ylabel('Amplitude'); grid on;

figure(2);
subplot(2,1,1);
stem(txBits(1:min(100, length(txBits))), 'filled');
title('Transmitted Bits'); ylim([-0.5 1.5]); grid on;
subplot(2,1,2);
stem(rxBits(1:min(100, length(rxBits))), 'filled');
title('Received Bits'); ylim([-0.5 1.5]); grid on;

[I_grid, Q_grid] = meshgrid(I_levels, Q_levels);
idealPoints = (I_grid(:) + 1j * Q_grid(:));

figure(3);
scatter(real(idealPoints), imag(idealPoints), 50, 'bo', 'filled'); hold on;
scatter(real(rx_symbols)*sqrt(42), imag(rx_symbols)*sqrt(42), 150, 'rx');
title('64QAM Constellation Diagram');
xlabel('In-phase (I)'); ylabel('Quadrature (Q)');
xline(0, 'k'); yline(0, 'k');
xlim([-8 8]); ylim([-8 8]); axis square; grid on;
