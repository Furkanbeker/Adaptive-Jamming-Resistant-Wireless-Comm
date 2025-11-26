%kanal bozulması yok

function LDPC_and_CRC_sim()

clc; close all;

% 1. input.txt'den veriyi oku
txt = fileread('input.txt');
if isempty(txt)
    error('input.txt dosyası boş. Lütfen içine bir metin girin.');
end
txBits = reshape(dec2bin(txt,8).' - '0', [], 1);  % sütun vektörü

% 2. CRC-32 ekle
crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
txBitsCRC = crcGen(txBits);

% 3. LDPC konfigürasyonu (DVB-S.2, oran 1/2)
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);

K = encCfg.NumInformationBits;
N = encCfg.BlockLength;

% 4. Girdi verisini K uzunluğa pad et
pad = mod(K - mod(length(txBitsCRC), K), K);
if pad ~= 0
    txBitsCRC = [txBitsCRC; zeros(pad, 1)];
end

% 5. LDPC encode (blok blok)
numBlocks = length(txBitsCRC) / K;
encodedBits = zeros(numBlocks * N, 1);
idx = 1;
for i = 1:K:length(txBitsCRC)
    block = txBitsCRC(i:i+K-1);
    encodedBlock = ldpcEncode(block, encCfg);
    encodedBits(idx:idx+N-1) = encodedBlock;
    idx = idx + N;
end

% 6. (Opsiyonel) Kanal bozulması
receivedBits = encodedBits;
% receivedBits(5000) = ~receivedBits(5000);  % test için bozma eklersen aç

% 7. LDPC decode (LLR ile, blok blok)
decodedBits = zeros(numBlocks * K, 1);
idx = 1;
for i = 1:N:length(receivedBits)
    codeword = receivedBits(i:i+N-1);
    llrInput = 1 - 2 * double(codeword);  % 0 -> +1, 1 -> -1
    decodedBlock = ldpcDecode(llrInput, decCfg, 50);
    decodedBits(idx:idx+K-1) = decodedBlock;
    idx = idx + K;
end

% 8. CRC kontrol
crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');
[rxBits, errFlag] = crcDet(decodedBits);

if errFlag
    error('❌ CRC kontrolü başarısız. Veri reddedildi.');
else
    disp('✅ CRC kontrolü başarılı. Veri doğru alındı.');
end

% 9. Padding varsa çıkar ve karaktere çevir
if pad ~= 0
    rxBits = rxBits(1:end-pad);
end
rxBitsMatrix = reshape(rxBits, 8, []).';
bitStrings = join(string(rxBitsMatrix), '', 2);
rxChars = char(bin2dec(bitStrings))';

% 10. output.txt dosyasına yaz
fid_out = fopen('output.txt', 'w');
fwrite(fid_out, rxChars, 'char');
fclose(fid_out);

disp("📄 output.txt dosyasına veri başarıyla yazıldı.");

end
