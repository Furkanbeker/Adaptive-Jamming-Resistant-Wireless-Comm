clc; close all; clear all;

SNR_dB_range = 0:2:20;
numSNR = length(SNR_dB_range);
numBits = 7200;

ber_bpsk = zeros(1, numSNR);
ber_qpsk = zeros(1, numSNR);
ber_16qam = zeros(1, numSNR);
ber_64qam = zeros(1, numSNR);

% BPSK
data = randi([0 1], 1, numBits);
symbols = 2 * data - 1;

for i = 1:numSNR
    rx = awgn(symbols, SNR_dB_range(i), 'measured');
    rxBits = double(rx > 0);
    ber_bpsk(i) = sum(rxBits ~= data) / numBits;
end

% QPSK (4-QAM)
data = randi([0 1], 1, numBits);
data_NRZ = 2 * data - 1;
s_p_data = reshape(data_NRZ, 2, []);
symbols = (s_p_data(1,:) + 1i * s_p_data(2,:)) / sqrt(2);

for i = 1:numSNR
    rx = awgn(symbols, SNR_dB_range(i), 'measured');
    rx_I = real(rx) > 0;
    rx_Q = imag(rx) > 0;
    rxBits = zeros(1, numBits);
    rxBits(1:2:end) = rx_I;
    rxBits(2:2:end) = rx_Q;
    ber_qpsk(i) = sum(rxBits ~= data) / numBits;
end

% 16-QAM
data = randi([0 1], 1, numBits);
bitGroups = reshape(data, 4, []).';
x = [-3 -1 3 1]; y = [3 1 -3 -1];
[X, Y] = meshgrid(x, y);
symMap = (X(:) + 1i*Y(:)) / sqrt(10);
decVals = bi2de(bitGroups, 'left-msb');
symbols = symMap(decVals + 1).';

for i = 1:numSNR
    rx = awgn(symbols, SNR_dB_range(i), 'measured');
    rxRep = repmat(rx.', 1, length(symMap));
    symRep = repmat(symMap.', length(rx), 1);
    dists = abs(rxRep - symRep);
    [~, minIdx] = min(dists, [], 2);
    rxBits = de2bi(minIdx - 1, 4, 'left-msb')';
    rxBits = rxBits(:).';
    ber_16qam(i) = sum(rxBits ~= data) / numBits;
end

% 64-QAM
data = randi([0 1], 1, numBits);
bitGroups = reshape(data, 6, []).';
i_gray = [-7 -5 -1 -3 7 5 1 3];
q_gray = [7 5 1 3 -7 -5 -1 -3];
[X, Y] = meshgrid(i_gray, q_gray);
symMap = (X(:) + 1i*Y(:)) / sqrt(42);
decVals = bi2de(bitGroups, 'left-msb');
symbols = symMap(decVals + 1).';

for i = 1:numSNR
    rx = awgn(symbols, SNR_dB_range(i), 'measured');
    rxRep = repmat(rx.', 1, length(symMap));
    symRep = repmat(symMap.', length(rx), 1);
    dists = abs(rxRep - symRep);
    [~, minIdx] = min(dists, [], 2);
    rxBits = de2bi(minIdx - 1, 6, 'left-msb')';
    rxBits = rxBits(:).';
    ber_64qam(i) = sum(rxBits ~= data) / numBits;
end

% Results
fprintf('\n--- BER Comparison ---\n');
fprintf('%6s | %10s | %10s | %10s | %10s\n', 'SNR(dB)', 'BPSK', 'QPSK', '16QAM', '64QAM');
fprintf(repmat('-',1,60)); fprintf('\n');
for i = 1:numSNR
    fprintf('%6d | %10.5f | %10.5f | %10.5f | %10.5f\n', ...
        SNR_dB_range(i), ber_bpsk(i), ber_qpsk(i), ber_16qam(i), ber_64qam(i));
end

% Plot
figure;
semilogy(SNR_dB_range, ber_bpsk, '-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB_range, ber_qpsk, '-s', 'LineWidth', 2);
semilogy(SNR_dB_range, ber_16qam, '-^', 'LineWidth', 2);
semilogy(SNR_dB_range, ber_64qam, '-d', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title('BER Comparison of Modulation Schemes');
legend('BPSK', 'QPSK', '16-QAM', '64-QAM', 'Location', 'southwest');
