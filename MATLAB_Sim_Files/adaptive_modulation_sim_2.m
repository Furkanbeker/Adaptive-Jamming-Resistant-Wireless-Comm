function adaptive_modulation()

clc; close all;

% SNR değerini al
EbNo_dB = input('SNR (dB) değerini girin: ');

% SNR değerine göre modülasyon seçimi
if EbNo_dB < 4
    modType = 'BPSK';
elseif EbNo_dB < 8
    modType = '4QAM';
elseif EbNo_dB < 15
    modType = '16QAM';
else
    modType = '64QAM';
end

fprintf('Seçilen modülasyon tipi: %s\n', modType);

% Giriş verisini oku
txt = fileread('input.txt');
txBits = reshape(dec2bin(txt, 8).' - '0', 1, []);
L = getBitsPerSymbol(modType);

% Padding
pad = mod(L - mod(length(txBits), L), L);
if pad ~= 0
    txBits = [txBits, zeros(1, pad)];
end

% Modülasyon
[symbols, normFactor] = modulate_bits(txBits, modType);

% Gürültü ekle
Es = mean(abs(symbols).^2);
Eb = Es / L;
N0 = Eb / (10^(EbNo_dB/10));
noise = sqrt(N0/2) * (randn(1,length(symbols)) + 1j*randn(1,length(symbols)));
rx_symbols = symbols + noise;

% Demodülasyon
rxBits = demodulate_symbols(rx_symbols, modType, normFactor);

% Uzunluk eşitleme (byte'a bölünecek kadar)
minLen = floor(min(length(txBits), length(rxBits)) / 8) * 8;
txBits = txBits(1:minLen);
rxBits = rxBits(1:minLen);

% BER hesapla
ber = sum(rxBits ~= txBits) / minLen;
fprintf('Eb/No = %.1f dB için Bit Error Rate (BER): %.5f\n', EbNo_dB, ber);

% Karaktere dönüştür ve dosyaya yaz
rxBitsMatrix = reshape(double(rxBits), 8, []).';   % <-- Burada double dönüşüm ekliyoruz
bitStrings = join(string(rxBitsMatrix), '', 2);    % 8-bit string üret
rxChars = char(bin2dec(bitStrings))';              % ASCII karaktere dönüştür
fid_out = fopen('output.txt', 'w');
fwrite(fid_out, rxChars, 'char');
fclose(fid_out);

% Bit karşılaştırma grafiği
figure;
subplot(2,1,1);
stem(txBits(1:min(100, end)), 'filled'); title('Gönderilen Bitler'); ylim([-0.5 1.5]); grid on;
subplot(2,1,2);
stem(rxBits(1:min(100, end)), 'filled'); title('Alınan Bitler'); ylim([-0.5 1.5]); grid on;

% Konstelasyon diyagramı
if strcmp(modType, 'BPSK')
    idealPoints = [-1 1];
    figure;
    scatter(real(idealPoints), zeros(1,length(idealPoints)), 50, 'bo', 'filled'); hold on;
    scatter(real(rx_symbols), zeros(1,length(rx_symbols)), 100, 'rx');
    title(['BPSK Constellation @ ', num2str(EbNo_dB), ' dB']);
    xlabel('In-Phase'); ylabel('Quadrature'); ylim([-0.5 0.5]);
    xline(0); yline(0); axis square; grid on;
else
    levels = get_levels(modType);
    idealPoints = [];
    for i = 1:length(levels)
        for j = 1:length(levels)
            idealPoints(end+1) = (levels(i) + 1j * levels(j)) / normFactor;
        end
    end
    figure;
    scatter(real(idealPoints), imag(idealPoints), 50, 'bo', 'filled'); hold on;
    scatter(real(rx_symbols), imag(rx_symbols), 100, 'rx');
    title(['Constellation: ', modType, ' @ ', num2str(EbNo_dB), ' dB']);
    xlabel('In-Phase'); ylabel('Quadrature');
    xline(0); yline(0); axis square; grid on;
end

end

% Bit başına sembol
function L = getBitsPerSymbol(modType)
switch modType
    case 'BPSK'; L = 1;
    case '4QAM'; L = 2;
    case '16QAM'; L = 4;
    case '64QAM'; L = 6;
    otherwise; error('Unknown modulation type');
end
end

% Seviye dizileri
function levels = get_levels(modType)
switch modType
    case 'BPSK'; levels = [-1 1];
    case '4QAM'; levels = [-1 1];
    case '16QAM'; levels = [-3 -1 1 3];
    case '64QAM'; levels = [-7 -5 -3 -1 1 3 5 7];
end
end

% Bitleri modüle et
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

% Sembolden bit çözümleri
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
