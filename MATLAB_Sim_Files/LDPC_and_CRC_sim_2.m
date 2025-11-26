function LDPC_CRC_AWGN_Full_Compare()

clc; close all;

% 1. input.txt'den veri oku
txt = fileread('input.txt');
if isempty(txt)
    error('input.txt dosyası boş. Lütfen içine bir metin girin.');
end
txBits = reshape(dec2bin(txt,8).' - '0', [], 1);

% 2. CRC-32 ekle
crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
txBitsCRC = crcGen(txBits);

% 3. LDPC konfigürasyonu
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
K = encCfg.NumInformationBits;
N = encCfg.BlockLength;

% 4. Padding işlemi
pad = mod(K - mod(length(txBitsCRC), K), K);
if pad ~= 0
    txBitsCRC = [txBitsCRC; zeros(pad, 1)];
end

% 5. LDPC encode
numBlocks = length(txBitsCRC) / K;
encodedBits = zeros(numBlocks * N, 1);
idx = 1;
for i = 1:K:length(txBitsCRC)
    block = txBitsCRC(i:i+K-1);
    encodedBlock = ldpcEncode(block, encCfg);
    encodedBits(idx:idx+N-1) = encodedBlock;
    idx = idx + N;
end

% 6. BPSK modülasyon
bpskSymbols = 1 - 2 * encodedBits;

% 7. AWGN kanal
EbNo_dB = 1.0;  % SNR
L = 1;
Es = mean(abs(bpskSymbols).^2);
Eb = Es / L;
N0 = Eb / (10^(EbNo_dB/10));
noise = sqrt(N0/2) * randn(length(bpskSymbols),1);
rxSymbols = bpskSymbols + noise;

% 8. LLR hesapla
llrInput = (2 * rxSymbols) / N0;

% 9. LDPC decode
decodedBits = zeros(numBlocks * K, 1);
idx = 1;
for i = 1:N:length(llrInput)
    llrBlock = llrInput(i:i+N-1);
    decodedBlock = ldpcDecode(llrBlock, decCfg, 50);
    decodedBits(idx:idx+K-1) = decodedBlock;
    idx = idx + K;
end

% 10. CRC kontrol
crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');
[rxBits, errFlag] = crcDet(decodedBits);

% 11. BER hesaplaması
hardDecoded = double(rxSymbols < 0);  % LDPC öncesi bit tahmini
preLDPCBits = hardDecoded(1:length(txBitsCRC));
ber_preLDPC = mean(preLDPCBits ~= txBitsCRC);
ber_postLDPC = mean(decodedBits(1:length(txBitsCRC)) ~= txBitsCRC);

fprintf('🔧 Bozulan bit sayısı (LDPC öncesi): %d\n', sum(preLDPCBits ~= txBitsCRC));
fprintf('📉 BER (LDPC öncesi): %.4f\n', ber_preLDPC);
fprintf('📈 BER (LDPC sonrası): %.4f\n', ber_postLDPC);

% 12. CRC sonucu
if errFlag
    disp('❌ CRC kontrolü başarısız. Veri reddedildi.');
else
    disp('✅ CRC kontrolü başarılı. Veri doğru alındı.');
end

% 13. output_1.txt → LDPC öncesi (hard decision)
if pad ~= 0
    preLDPCBits = preLDPCBits(1:end-pad);
end
usableLength = floor(length(preLDPCBits)/8)*8;
preLDPCBits = preLDPCBits(1:usableLength);
preMatrix = reshape(preLDPCBits, 8, []).';

preChars = char(bin2dec(join(string(preMatrix), '', 2)))';
fid1 = fopen('output_1.txt', 'w');
fwrite(fid1, preChars, 'char');
fclose(fid1);

% 14. output_2.txt → LDPC sonrası (decoded)
if pad ~= 0
    rxBits = rxBits(1:end-pad);
end
rxMatrix = reshape(rxBits, 8, []).';
rxChars = char(bin2dec(join(string(rxMatrix), '', 2)))';
fid2 = fopen('output_2.txt', 'w');
fwrite(fid2, rxChars, 'char');
fclose(fid2);

disp("📁 output_1.txt → LDPC öncesi (bozulmuş)");
disp("📁 output_2.txt → LDPC sonrası (düzeltilmiş)");

% -------------------- GRAFİKLER --------------------

% A) BPSK Sinyal + Gürültü
figure;
plot(bpskSymbols, 'b'); hold on;
plot(rxSymbols, 'r'); title('BPSK Sinyali ve AWGN Sonrası');
legend('Gönderilen','Gelen'); grid on;

% B) LLR Histogram
figure;
histogram(llrInput, 100);
title('LLR Değer Dağılımı'); xlabel('LLR'); ylabel('Sayı'); grid on;

% C) İlk 100 bit karşılaştırması
figure;
stem(txBits(1:100), 'b'); hold on;
stem(rxBits(1:100), 'r'); title('Gönderilen vs Alınan Bitler (İlk 100)');
legend('Tx','Rx'); ylim([-0.5 1.5]); grid on;

end
