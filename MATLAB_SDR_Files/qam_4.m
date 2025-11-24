clc; clear; close all;

% File Paths
inputPath = 'C:\Users\mfb36\Desktop\input.txt';
outputPath = 'C:\Users\mfb36\Desktop\output.txt';

% Communication Parameters
br = 5e4;
fs = 1e6;
fc = 915e6;
samplesPerSymbol = fs / br;
txID = 'usb:1.0';
rxID = 'usb:1.1';

% Read Input Text
fid = fopen(inputPath, 'r');
if fid == -1
    error('input.txt could not be found!');
end
textToSend = fread(fid, '*char')';
fclose(fid);

% Text to QAM Symbols
binData = dec2bin(textToSend, 8);
txBits = reshape(binData.' - '0', 1, []);
bitPairs = reshape(txBits, 2, []).';
decVals = bi2de(bitPairs, 'left-msb');
qamSymbols = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);
modSignal = qamSymbols(decVals + 1);

% Preamble
preambleBits = [1 0 1 0 0 1 0 1];
preamblePairs = reshape(preambleBits, 2, []).';
preambleVals = bi2de(preamblePairs, 'left-msb');
preambleSymbols = qamSymbols(preambleVals + 1);

% Upsample
txSymbols = [preambleSymbols modSignal];
txUpsampled = repelem(txSymbols, samplesPerSymbol);

% PlutoSDR Transmitter
tx = sdrtx('Pluto', ...
    'CenterFrequency', fc, ...
    'BasebandSampleRate', fs, ...
    'Gain', 0, ...
    'RadioID', txID);

release(tx);
tx.transmitRepeat(txUpsampled.');

pause(2);

% PlutoSDR Receiver
rx = sdrrx('Pluto', ...
    'CenterFrequency', fc, ...
    'BasebandSampleRate', fs, ...
    'SamplesPerFrame', 5e5, ...
    'Gain', 35, ...
    'OutputDataType', 'double', ...
    'RadioID', rxID);

disp("Receiving...");
rxData = rx();
rxData = rxData(:);

% Add AWGN Noise
addNoise = true;
noisePowerInput = 0.2;
if addNoise
    noise = sqrt(noisePowerInput/2) * (randn(size(rxData)) + 1j * randn(size(rxData)));
    rxData = rxData + noise;
    fprintf("Manual AWGN added. Input power = %.6f\n", noisePowerInput);
else
    noise = zeros(size(rxData));
end

% Downsample
rxDownsampled = buffer(rxData, samplesPerSymbol);
rxSymbols = mean(rxDownsampled, 1);

% Synchronization
corr = abs(conv(rxSymbols, fliplr(conj(preambleSymbols)), 'valid'));
[~, startIdx] = max(corr);
alignedStart = startIdx + length(preambleSymbols);

payloadLen = length(modSignal);
if alignedStart + payloadLen - 1 > length(rxSymbols)
    disp(['rxSymbols length: ' num2str(length(rxSymbols))]);
    disp(['alignedStart: ' num2str(alignedStart)]);
    disp(['payloadLen: ' num2str(payloadLen)]);
    error('Synchronization failed: Not enough symbols after preamble. Try increasing SamplesPerFrame or check your signal.');
end

% Extract Payload
rxPayload = rxSymbols(alignedStart : alignedStart + payloadLen - 1);

% Amplitude Normalization
ampEst = mean(abs(rxPayload)) / mean(abs(modSignal));
rxPayload_norm = rxPayload / ampEst;

% Carrier Phase Recovery
angleGrid = linspace(0, 2*pi, 128);
minErr = inf;
for theta = angleGrid
    rotated = rxPayload_norm * exp(-1j * theta);
    rxBitsTmp = reshape(de2bi(qamDemodulate(rotated), 2, 'left-msb').', 1, []);
    err = sum(rxBitsTmp ~= txBits);
    if err < minErr
        minErr = err;
        bestAngle = theta;
        rxBitsBest = rxBitsTmp;
        rxPayloadAligned = rotated;
    end
end

% Bit and Char Decoding
rxBits = rxBitsBest;
rxChar = char(bin2dec(num2str(reshape(rxBits, 8, []).')));

% BER Calculation
ber = sum(rxBits ~= txBits) / length(txBits);

% Write Decoded Output to File
fid = fopen(outputPath, 'w');
if fid == -1
    error('output.txt could not be written!');
end
fprintf(fid, '%s', rxChar);
fclose(fid);
fprintf("Decoded output saved to: %s\n", outputPath);

% SNR, EVM and Drift Calculation
noiseEstimate = rxPayloadAligned - modSignal;
signalPower = mean(abs(modSignal).^2);
noisePower = mean(abs(noiseEstimate).^2);
snr = 10 * log10(signalPower / noisePower);
evmRMS = sqrt(mean(abs(noiseEstimate).^2));
evmPercent = evmRMS * 100;
drift = mean(abs(rxPayloadAligned));

fprintf('\nTransmitted: %s\n', textToSend);
fprintf('Received   : %s\n', rxChar);
fprintf('BER        : %.4f\n', ber);
fprintf('SNR        : %.2f dB\n', snr);
fprintf('EVM        : %.2f %% RMS\n', evmPercent);
fprintf('Drift      : %.3f\n', drift);

% Constellation Diagram
figure;
hold on;
plot(real(qamSymbols), imag(qamSymbols), 'ko', 'MarkerSize', 16, 'LineWidth', 2);
scatter(real(modSignal), imag(modSignal), 120, 'bo', 'LineWidth', 2);
scatter(real(rxPayloadAligned), imag(rxPayloadAligned), 120, 'rx', 'LineWidth', 2);
legend('4QAM Reference', 'Transmitted', 'Received');
title('Constellation Diagram');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-1.5 1.5 -1.5 1.5]); grid on;

% Bitstream Comparison Plot
N = min(length(rxBits), length(txBits));
figure;
subplot(2,1,1); stairs(txBits(1:N), 'b', 'LineWidth', 1.2); title('Transmitted Bitstream'); ylim([-0.2 1.2]);
subplot(2,1,2); stairs(rxBits(1:N), 'r', 'LineWidth', 1.2); title('Received Bitstream'); ylim([-0.2 1.2]); xlabel('Bit Index');

% Power Spectral Density Plot
figure;
pwelch(rxData, [], [], [], fs, 'centered');
title('Received Signal Power Spectral Density');

% QPSK/4QAM Demodulation Function
function out = qamDemodulate(symbols)
    ref = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);
    out = zeros(1, length(symbols));
    for k = 1:length(symbols)
        [~, idx] = min(abs(symbols(k) - ref));
        out(k) = idx - 1;
    end
end
