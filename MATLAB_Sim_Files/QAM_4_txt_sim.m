clc; clear all; close all;

% File Path Settings
inputPath = 'C:\Users\mfb36\Desktop\input.txt';
outputPath = 'C:\Users\mfb36\Desktop\output.txt';

% Read text from input.txt
fid = fopen(inputPath, 'r');
if fid == -1
    error('Unable to open input.txt. Check the file path.');
end
txt = fread(fid, '*char')';
fclose(fid);

% ASCII to Bit Sequence Conversion
decData = dec2bin(txt, 8);
txBits = reshape(decData.' - '0', 1, []);
numBits = length(txBits);
if mod(numBits,2) ~= 0
    txBits(end+1) = 0;
    numBits = numBits + 1;
end

% Modulation (4QAM)
symbols = zeros(1, numBits/2);
for i = 1:2:numBits
    b = txBits(i:i+1);
    if isequal(b, [0 0])
        s = 1 + 1j;
    elseif isequal(b, [0 1])
        s = -1 + 1j;
    elseif isequal(b, [1 1])
        s = -1 - 1j;
    else 
        s = 1 - 1j;
    end
    symbols((i+1)/2) = s;
end
symbols = symbols / sqrt(2);  % Normalize power

% Noise addition and SNR calculation
signal_power = mean(abs(symbols).^2);
noise_power = 0.1; 
snr_linear = signal_power / noise_power;
SNR_dB = 10 * log10(snr_linear);

noise_sym = sqrt(noise_power/2) * (randn(1,length(symbols)) + 1j*randn(1,length(symbols)));
rx_symbols = symbols + noise_sym;

% Demodulation
rxBits = zeros(1, numBits);
for i = 1:length(rx_symbols)
    I = real(rx_symbols(i));
    Q = imag(rx_symbols(i));

    if I > 0 && Q > 0
        rxBits(2*i-1:2*i) = [0 0];
    elseif I < 0 && Q > 0
        rxBits(2*i-1:2*i) = [0 1];
    elseif I < 0 && Q < 0
        rxBits(2*i-1:2*i) = [1 1];
    else
        rxBits(2*i-1:2*i) = [1 0];
    end
end

% BER Calculation
excess = mod(length(rxBits), 8);
if excess ~= 0
    rxBits = rxBits(1:end-excess);
    txBits = txBits(1:end-excess);
end
ber = sum(rxBits ~= txBits) / length(rxBits);

% Display metrics
fprintf('Signal Power: %.4f\n', signal_power);
fprintf('Noise Power: %.4f\n', noise_power);
fprintf('Calculated SNR: %.2f dB\n', SNR_dB);
fprintf('Bit Error Rate (BER): %.4f\n', ber);

% Write received characters to output file
rxBitsMatrix = reshape(rxBits, 8, []).';
rxChars = char(bin2dec(num2str(rxBitsMatrix)))';
fid_out = fopen(outputPath, 'w');
fwrite(fid_out, rxChars);
fclose(fid_out);

% Physical waveform generation (optional visualization)
br = 1e6; f = br; T = 1/br;
fs = 100 * f;
t = 0:1/fs:T-1/fs;
samplesPerSymbol = length(t);
y_mod = []; y_in = []; y_qd = [];

for i = 1:length(symbols)
    I = real(symbols(i));
    Q = imag(symbols(i));
    in = I * cos(2*pi*f*t);
    qd = Q * sin(2*pi*f*t);
    y_in = [y_in in];
    y_qd = [y_qd qd];
    y_mod = [y_mod in + qd];
end

% Add time-domain noise
time_noise = sqrt(noise_power) * randn(1, length(y_mod));
y_noisy = y_mod + time_noise;
tt = linspace(0, T*length(symbols), length(y_mod));

% Visualization
figure(1);
subplot(4,1,1);
plot(tt, y_in, 'LineWidth', 1.2);
title('4QAM In-phase Component (I)');
ylabel('Amplitude'); grid on;

subplot(4,1,2);
plot(tt, y_qd, 'LineWidth', 1.2);
title('4QAM Quadrature Component (Q)');
ylabel('Amplitude'); grid on;

subplot(4,1,3);
plot(tt, time_noise, 'LineWidth', 1.2);
title('Channel Noise');
ylabel('Amplitude'); grid on;

subplot(4,1,4);
plot(tt, y_noisy, 'r', 'LineWidth', 1.2);
title('Noisy 4QAM Signal');
xlabel('Time'); ylabel('Amplitude'); grid on;

figure(2);
subplot(2,1,1); stem(txBits(1:min(100, length(txBits))), 'filled'); title('Transmitted Bits'); ylim([-0.5 1.5]); grid on;
subplot(2,1,2); stem(rxBits(1:min(100, length(rxBits))), 'filled'); title('Received Bits'); ylim([-0.5 1.5]); grid on;

figure(3);
scatter(real(symbols), imag(symbols), 50, 'bo', 'filled'); hold on;
scatter(real(rx_symbols), imag(rx_symbols), 150, 'rx');
title('4QAM Constellation Diagram');
xlabel('In-phase (I)'); ylabel('Quadrature (Q)');
xline(0, 'k'); yline(0, 'k');
xlim([-2 2]); ylim([-2 2]); axis square; grid on;