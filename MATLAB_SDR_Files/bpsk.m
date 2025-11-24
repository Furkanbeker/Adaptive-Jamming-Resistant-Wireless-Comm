clc; clear; close all;

inputPath = 'C:\Users\mfb36\Desktop\input.txt';
outputPath = 'C:\Users\mfb36\Desktop\output.txt';

% Parameters
br = 1e5;
fs = 10e6;
fc = 915e6;
T = 1 / br;
t = 0:1/fs:T-1/fs;
samplesPerSymbol = length(t);

fid = fopen(inputPath, 'r');
if fid == -1
    error('input.txt could not find!');
end
textToSend = fread(fid, '*char')';
fclose(fid);

% Convert text to binary sequence
binData = dec2bin(char(textToSend), 8);
txBits = reshape(binData.' - '0', 1, []);
symbols = 2 * txBits - 1;

% Define preamble
preambleBits = [1 0 1 1 0 0 1 0];
preambleSymbols = 2 * preambleBits - 1;

% BPSK Modulation (with carrier)
modulatedSignal = [];
for i = 1:length(preambleSymbols)
    modulatedSignal = [modulatedSignal preambleSymbols(i)*cos(2*pi*fc*t)];
end
for i = 1:length(symbols)
    modulatedSignal = [modulatedSignal symbols(i)*cos(2*pi*fc*t)];
end
txSignal = [zeros(1, samplesPerSymbol*100), modulatedSignal];
txSignal = txSignal / max(abs(txSignal));

% Pluto Transmitter
tx = sdrtx('Pluto', ...
    'CenterFrequency', fc, ...
    'BasebandSampleRate', fs, ...
    'Gain', -10, ...
    'RadioID', 'usb:1.3.5');
release(tx);
tx.transmitRepeat(complex(txSignal.'));

% Pluto Receiver
rx = sdrrx('Pluto', ...
    'CenterFrequency', fc, ...
    'BasebandSampleRate', fs, ...
    'SamplesPerFrame', 2e5, ...
    'OutputDataType','double', ...
    'RadioID', 'usb:1.40.5');

pause(1);
rxSignal = rx();
rxSignal = rxSignal(:);

% Add AWGN Noise
addNoise = true;
noisePowerInput = 0.5;
if addNoise
    noise = sqrt(noisePowerInput/2) * (randn(size(rxSignal)) + 1j * randn(size(rxSignal)));
    rxSignal = rxSignal + noise;
    fprintf("Manual AWGN added. Input power = %.6f\n", noisePowerInput);
else
    noise = zeros(size(rxSignal));
end

% Carrier removal (Demodulation)
n = (0:length(rxSignal)-1).';
rxBaseband = rxSignal .* cos(2*pi*fc*n/fs);
rxBaseband = lowpass(rxBaseband, br/2, fs);

% Preamble Correlation
basePreamble = [];
for i = 1:length(preambleSymbols)
    basePreamble = [basePreamble preambleSymbols(i)*ones(1, samplesPerSymbol)];
end
searchLength = min(3e5, length(rxBaseband));
corr = conv(rxBaseband(1:searchLength), fliplr(basePreamble), 'valid');
corr = abs(corr) / max(abs(corr));
[~, startIdx] = max(corr);
rxAligned = rxBaseband(startIdx + 1:end);

% Downsampling
totalBits = length(preambleBits) + length(txBits);
rxSamples = zeros(1, totalBits);
for i = 1:totalBits
    idx = (i-1)*samplesPerSymbol + 1;
    if idx + samplesPerSymbol - 1 <= length(rxAligned)
        rxSamples(i) = mean(rxAligned(idx : idx+samplesPerSymbol-1));
    else
        rxSamples(i) = 0;
    end
end

% Bit decision
rxData = rxSamples(length(preambleBits)+1 : end);
rxBits = rxData > 0;

% BER Calculation
rxBits = rxBits(1:length(txBits));
if sum(rxBits ~= txBits) > sum(~rxBits ~= txBits)
    rxBits = ~rxBits;
end
ber = sum(rxBits ~= txBits) / length(rxBits);

% ASCII Conversion
rxMatrix = reshape(rxBits, 8, []).';
rxChars = char(bin2dec(num2str(rxMatrix)))';

% Terminal Output and TXT
fprintf("\nTransmitted: %s\n", textToSend);
fprintf("Received:    %s\n", rxChars);
fprintf("BER:         %.4f\n", ber);

% Write to Output.txt
fid = fopen(outputPath, 'w');
if fid == -1
    error('output.txt could not be written!');
end
fprintf(fid, '%s', rxChars);
fclose(fid);
fprintf("Decoded output saved to: %s\n", outputPath);

% SNR Calculation (from known clean signal and added noise)
cleanSignal = 2*double(txBits) - 1;
signalPower = mean(cleanSignal.^2);                      
measuredNoisePower = mean(abs(noise).^2);                
snr = 10 * log10(signalPower / measuredNoisePower + eps);

fprintf("Average Signal Power: %.6f\n", signalPower);
fprintf("Measured Noise Power: %.6f\n", measuredNoisePower);
fprintf("SNR (from power ratio): %.2f dB\n", snr);

% Constellation Diagram (BPSK) - Only Crosses
txSymbols = symbols;
rxSymbols = 2*double(rxData > 0) - 1;
jitterAmount = 0.02;

figure; hold on;

correctIdx = find(rxBits == txBits);
errorIdx   = find(rxBits ~= txBits);

xCorr = rxSymbols(correctIdx) + randn(size(correctIdx)) * jitterAmount;
yCorr = randn(size(correctIdx)) * jitterAmount * 0.5;

xErr = rxSymbols(errorIdx) + randn(size(errorIdx)) * jitterAmount;
yErr = randn(size(errorIdx)) * jitterAmount * 0.5;

h1 = scatter(xCorr, yCorr, 100, 'g', 'x', 'LineWidth', 2, 'DisplayName', 'Correct Bits');
h2 = scatter(xErr, yErr, 120, 'r', 'x', 'LineWidth', 2, 'DisplayName', 'Bit Errors');

legend([h1 h2], 'Location', 'northeast');
title('Constellation Diagram (BPSK) - Cross Only');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2 2 -1 1]); grid on; box on;

% Bitstream Comparison Plot
N = min(length(rxBits), length(txBits));
figure;
subplot(2,1,1); stairs(txBits(1:N), 'b', 'LineWidth', 1.2); title('Transmitted Bitstream'); ylim([-0.2 1.2]); grid on;
subplot(2,1,2); stairs(rxBits(1:N), 'r', 'LineWidth', 1.2); title('Received Bitstream'); ylim([-0.2 1.2]); xlabel('Bit Index'); grid on;

% Power Spectral Density
figure;
pwelch(rxSignal, [], [], [], fs, 'centered');
title('Received Signal Power Spectral Density');
