function adaptive_modulation_ldpc_crc_sim()
clc; close all;

% 1. SNR & modülasyon seçimi
EbNo_dB = input('SNR (dB) değerini girin: ');
modType = select_modulation(EbNo_dB);
fprintf('Seçilen modülasyon: %s\n', modType);

% 2. Girdi verisi
txt = fileread('input.txt');
if isempty(txt), error('input.txt boş!'); end
txBits = reshape(dec2bin(txt,8).' - '0', [], 1);

% 3. CRC ekle
crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
txCRC = crcGen(txBits);

% 4. LDPC ayarları
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
K = encCfg.NumInformationBits; N = encCfg.BlockLength;

% 5. Padding
pad = mod(K - mod(length(txCRC),K),K);
if pad > 0, txCRC = [txCRC; zeros(pad,1)]; end

% 6. LDPC encode
numBl = length(txCRC)/K;
encBits = zeros(numBl*N,1);
for i=1:numBl
    encBits((i-1)*N+1:i*N) = ldpcEncode(txCRC((i-1)*K+1:i*K), encCfg);
end

% --- BPSK/QAM MODÜLASYON BLOĞU ---
switch modType
    case 'BPSK'
        L = 1; M = 2;
        % BPSK mapping ve bit tipi düzeltme
        if islogical(encBits), encBits = double(encBits); end
        symbols = 1 - 2*encBits;   % textbook BPSK: 0->+1, 1->-1
        Eb = 1;
        N0 = Eb / (10^(EbNo_dB/10));
        noise = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
        rxSym = symbols + noise;
        llr = (2/N0)*real(rxSym);
    case '4QAM'
        L = 2; M = 4;
        bits = reshape(encBits, L, []).';
        symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
        Es = mean(abs(symbols).^2); Eb = Es / L;
        N0 = Eb / (10^(EbNo_dB/10));
        noise = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
        rxSym = symbols + noise;
        llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
            'UnitAveragePower', true, 'NoiseVariance', N0);
    case '16QAM'
        L = 4; M = 16;
        bits = reshape(encBits, L, []).';
        symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
        Es = mean(abs(symbols).^2); Eb = Es / L;
        N0 = Eb / (10^(EbNo_dB/10));
        noise = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
        rxSym = symbols + noise;
        llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
            'UnitAveragePower', true, 'NoiseVariance', N0);
    case '64QAM'
        L = 6; M = 64;
        bits = reshape(encBits, L, []).';
        symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
        Es = mean(abs(symbols).^2); Eb = Es / L;
        N0 = Eb / (10^(EbNo_dB/10));
        noise = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
        rxSym = symbols + noise;
        llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
            'UnitAveragePower', true, 'NoiseVariance', N0);
    otherwise
        error('Bilinmeyen modülasyon!');
end

llr = llr(:);

% 8. BER LDPC öncesi
rxHard = llr < 0;
ber_pre = mean(rxHard ~= encBits);
fprintf('📉 BER (öncesi): %.5f\n', ber_pre);

% 9. LDPC decode
decBits = zeros(numBl*K,1);
for i=1:numBl
    decBits((i-1)*K+1:i*K) = ldpcDecode(llr((i-1)*N+1:i*N), decCfg, 50);
end

% 10. CRC kontrol
crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');
[rxCRC, err] = crcDet(decBits);
ber_post = mean(decBits(1:length(txCRC)) ~= txCRC);
fprintf('📈 BER (sonrası): %.5f\n', ber_post);

if err
    warning('❌ CRC başarısız!');
    rxCRC = decBits(1:end - pad); % fallback
else
    disp('✅ CRC başarılı.');
    if pad > 0 && length(rxCRC) >= pad
        rxCRC = rxCRC(1:end-pad);
    end
end

% 11. output_1.txt
write_bits_to_file('output_1.txt', rxHard);

% 12. output_2.txt
write_bits_to_file('output_2.txt', rxCRC);

% 13. Konstelasyon diyagramı
figure; 
scatter(real(symbols(1:min(1000,end))), imag(symbols(1:min(1000,end))), '.b'); hold on;
scatter(real(rxSym(1:min(1000,end))), imag(rxSym(1:min(1000,end))), '.r');
title(['Konstelasyon: ',modType,' @',num2str(EbNo_dB),' dB']); legend('Tx','Rx');

% 14. LLR histogramı
figure; histogram(llr,100); title('LLR Dağılımı');

% 15. İlk 100 bit karşılaştırması (altlı üstlü)
figure;
subplot(2,1,1); stem(txBits(1:100),'b','filled'); title('Gönderilen İlk 100 Bit'); ylim([-0.5 1.5]);
subplot(2,1,2); stem(rxCRC(1:min(100,end)),'r','filled'); title('Alınan İlk 100 Bit'); ylim([-0.5 1.5]);
end

function m = select_modulation(EbNo_dB)
    if EbNo_dB < 4, m = 'BPSK';
    elseif EbNo_dB < 8, m = '4QAM';
    elseif EbNo_dB < 15, m = '16QAM';
    else, m = '64QAM'; end
end

function write_bits_to_file(filename, bits)
    bits = bits(:);
    pad = mod(8 - mod(length(bits),8), 8);
    bits = [bits; zeros(pad,1)];
    bytes = reshape(bits, 8, []).';
    binStrs = num2str(bytes);
    binStrs = binStrs(:,~isspace(binStrs(1,:)));
    chars = char(bin2dec(binStrs))';
    fid = fopen(filename,'w');
    fwrite(fid, chars, 'char');
    fclose(fid);
end
