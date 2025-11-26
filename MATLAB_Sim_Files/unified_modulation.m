
function unified_modulation(modType, EbNo_dB)
% modType: 'BPSK', '4QAM', '16QAM', or '64QAM'
% EbNo_dB: SNR in dB

clc; close all;

% Input & Bit Generation
txt = fileread('C:\Users\mfb36\Desktop\input.txt');
txBits = reshape(dec2bin(txt, 8).' - '0', 1, []);
L = getBitsPerSymbol(modType);

% Pad bits
pad = mod(L - mod(length(txBits), L), L);
if pad ~= L
    txBits = [txBits, zeros(1, pad)];
end

% Modulation
[symbols, normFactor] = modulate_bits(txBits, modType);

% Noise
Es = mean(abs(symbols).^2);
Eb = Es / L;
N0 = Eb / (10^(EbNo_dB/10));
noise = sqrt(N0/2) * (randn(1,length(symbols)) + 1j*randn(1,length(symbols)));
rx_symbols = symbols + noise;

% Demodulation
rxBits = demodulate_symbols(rx_symbols, modType, normFactor);

% Remove padding
rxBits = rxBits(1:end - pad);
txBits = txBits(1:end - pad);

% BER Calculation
ber = sum(rxBits ~= txBits) / length(txBits);
fprintf('Modulation: %s\nEb/No: %.1f dB\nBER: %.5f\n', modType, EbNo_dB, ber);

% Output
rxBytes = char(bin2dec(reshape(char(rxBits+'0'), 8, []).')).';
fid = fopen('output.txt', 'w'); fwrite(fid, rxBytes); fclose(fid);

% Constellation Plot
scatterplot(rx_symbols);
title(['Constellation: ', modType, ' @ ', num2str(EbNo_dB), ' dB']);
grid on;

end

function L = getBitsPerSymbol(modType)
switch modType
    case 'BPSK'; L = 1;
    case '4QAM'; L = 2;
    case '16QAM'; L = 4;
    case '64QAM'; L = 6;
    otherwise; error('Unknown modulation type');
end
end

function [symbols, normFactor] = modulate_bits(bits, modType)
switch modType
    case 'BPSK'
        symbols = 2*bits - 1;
        normFactor = 1;
    case '4QAM'
        levels = [-1 1];
        normFactor = sqrt(2);
        symbols = [];
        for i = 1:2:length(bits)
            b = bits(i:i+1);
            I = levels(b(1)+1);
            Q = levels(b(2)+1);
            symbols(end+1) = (I + 1j*Q)/normFactor;
        end
    case '16QAM'
        levels = [-3 -1 1 3];
        normFactor = sqrt(10);
        symbols = [];
        for i = 1:4:length(bits)
            I = levels(bi2de(bits(i:i+1), 'left-msb')+1);
            Q = levels(bi2de(bits(i+2:i+3), 'left-msb')+1);
            symbols(end+1) = (I + 1j*Q)/normFactor;
        end
    case '64QAM'
        levels = [-7 -5 -3 -1 1 3 5 7];
        normFactor = sqrt(42);
        I_bin = de2bi(0:7, 3, 'left-msb');
        Q_bin = de2bi(0:7, 3, 'left-msb');
        symbols = [];
        for i = 1:6:length(bits)
            I_idx = bi2de(bits(i:i+2), 'left-msb')+1;
            Q_idx = bi2de(bits(i+3:i+5), 'left-msb')+1;
            I = levels(I_idx);
            Q = levels(Q_idx);
            symbols(end+1) = (I + 1j*Q)/normFactor;
        end
    otherwise
        error('Unsupported modulation');
end
end

function bits = demodulate_symbols(rx_symbols, modType, normFactor)
switch modType
    case 'BPSK'
        bits = real(rx_symbols) > 0;
    case '4QAM'
        bits = [];
        for r = rx_symbols
            I = real(r)*normFactor;
            Q = imag(r)*normFactor;
            bits = [bits, I>0, Q>0];
        end
    case '16QAM'
        levels = [-3 -1 1 3];
        bits = [];
        for r = rx_symbols
            I = real(r)*normFactor;
            Q = imag(r)*normFactor;
            [~, I_idx] = min(abs(levels - I));
            [~, Q_idx] = min(abs(levels - Q));
            bits = [bits, de2bi(I_idx-1,2,'left-msb'), de2bi(Q_idx-1,2,'left-msb')];
        end
    case '64QAM'
        levels = [-7 -5 -3 -1 1 3 5 7];
        bits = [];
        for r = rx_symbols
            I = real(r)*normFactor;
            Q = imag(r)*normFactor;
            [~, I_idx] = min(abs(levels - I));
            [~, Q_idx] = min(abs(levels - Q));
            bits = [bits, de2bi(I_idx-1,3,'left-msb'), de2bi(Q_idx-1,3,'left-msb')];
        end
    otherwise
        error('Unsupported modulation');
end
end