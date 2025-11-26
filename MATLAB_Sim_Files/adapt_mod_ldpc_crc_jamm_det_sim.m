function adapt_mod_ldpc_crc_jammer()
clc; close all;

%% 1. Kullanıcıdan SNR ve sistem parametreleri
EbNo_dB = input('SNR (dB) değerini girin: ');
center_freqs = [100e3 250e3 400e3];  % Kullanılabilir frekanslar
jammer_freq = 250e3;                 % Jammer burada!
jammer_amp = 2;                      % Jammer sinyal genliği
noise_level = 0.15;                  % Kanal gürültü genliği
fs = 1e6;                            % Örnekleme frekansı (Hz)
T = 1;                               % Sinyal süresi (sn)
t = 0:1/fs:T-1/fs;                   % Zaman vektörü

%% 2. Girdi verisi
txt = fileread('input.txt');
if isempty(txt), error('input.txt boş!'); end
txBits = reshape(dec2bin(txt,8).' - '0', [], 1);

%% 3. CRC ekle
crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
txCRC = crcGen(txBits);

%% 4. LDPC ayarları
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
K = encCfg.NumInformationBits; N = encCfg.BlockLength;

%% 5. Padding
pad = mod(K - mod(length(txCRC),K),K);
if pad > 0, txCRC = [txCRC; zeros(pad,1)]; end

%% 6. LDPC encode
numBl = length(txCRC)/K;
encBits = zeros(numBl*N,1);
for i=1:numBl
    encBits((i-1)*N+1:i*N) = ldpcEncode(txCRC((i-1)*K+1:i*K), encCfg);
end

%% 7. Frekans seçimi & JAMMER DETECTION + FREKANS ATLAMA
fprintf('Frekans taraması başlıyor...\n');
clean_centers = [];
for cf = center_freqs
    if abs(cf - jammer_freq) < 1
        jammer = jammer_amp*sin(2*pi*jammer_freq*t);
    else
        jammer = zeros(size(t));
    end
    noise = noise_level*randn(size(t));
    rx_total = jammer + noise;   % SİNYAL YOK!

    % --- GELİŞMİŞ JAMMER DETECTION ---
    Nfft = 2^nextpow2(length(rx_total));
    f = fs*(0:(Nfft/2))/Nfft;
    RX = fft(rx_total, Nfft);
    PSD = abs(RX/Nfft).^2;
    PSD = PSD(1:Nfft/2+1);

    margin = 500; % 500 Hz pencere
    idx = find(abs(f - cf) < margin);

    % En yüksek 5 değerin ortalamasını al
    [~, sortidx] = sort(PSD(idx), 'descend');
    top_vals = PSD(idx(sortidx(1:min(5,end))));
    power_dB = 10*log10(mean(top_vals));

    % Arka plan: bandın ilk %25'inin ortalaması
    bg_band = PSD(1:round(length(PSD)*0.25));
    bg_dB = 10*log10(mean(bg_band));
    threshold = bg_dB + 15;  % 15 dB üstü jammed

    fprintf('%.0f kHz: Güç = %.2f dB, BG = %.2f dB, Eşik = %.2f dB -> ', cf/1e3, power_dB, bg_dB, threshold);
    if power_dB > threshold
        fprintf('JAMMED\n');
    else
        fprintf('TEMİZ\n');
        clean_centers = [clean_centers, cf];
    end
end

%% 8. Frekans atlama kararı
if isempty(clean_centers)
    error('Tüm merkez frekanslar jammed! Kaçacak yer yok.');
else
    new_center = clean_centers(1); % Temiz ilk frekansı seç
    fprintf('\nKullanılacak frekans: %.0f kHz\n', new_center/1e3);
end

%% 9. Yeni merkezde iletim (temiz frekans)
signal = sin(2*pi*new_center*t);
if abs(new_center - jammer_freq) < 1
    jammer = jammer_amp*sin(2*pi*jammer_freq*t);
else
    jammer = zeros(size(t));
end
noise = noise_level*randn(size(t));
rx_total = signal + jammer + noise;

% --- Modülasyon seçimi (SNR'a göre) ---
modType = select_modulation(EbNo_dB);
fprintf('Seçilen modülasyon: %s\n', modType);

% --- MODÜLASYON BLOĞU ---
switch modType
    case 'BPSK'
        L = 1; M = 2;
        if islogical(encBits), encBits = double(encBits); end
        symbols = 1 - 2*encBits;
        Eb = 1;
        N0 = Eb / (10^(EbNo_dB/10));
        noise_ch = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
        rxSym = symbols + noise_ch;
        llr = (2/N0)*real(rxSym);
    case '4QAM'
        L = 2; M = 4;
        bits = reshape(encBits, L, []).';
        symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
        Es = mean(abs(symbols).^2); Eb = Es / L;
        N0 = Eb / (10^(EbNo_dB/10));
        noise_ch = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
        rxSym = symbols + noise_ch;
        llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
            'UnitAveragePower', true, 'NoiseVariance', N0);
    case '16QAM'
        L = 4; M = 16;
        bits = reshape(encBits, L, []).';
        symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
        Es = mean(abs(symbols).^2); Eb = Es / L;
        N0 = Eb / (10^(EbNo_dB/10));
        noise_ch = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
        rxSym = symbols + noise_ch;
        llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
            'UnitAveragePower', true, 'NoiseVariance', N0);
    case '64QAM'
        L = 6; M = 64;
        bits = reshape(encBits, L, []).';
        symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
        Es = mean(abs(symbols).^2); Eb = Es / L;
        N0 = Eb / (10^(EbNo_dB/10));
        noise_ch = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
        rxSym = symbols + noise_ch;
        llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
            'UnitAveragePower', true, 'NoiseVariance', N0);
    otherwise
        error('Bilinmeyen modülasyon!');
end
llr = llr(:);

%% 10. BER LDPC öncesi
rxHard = llr < 0;
ber_pre = mean(rxHard ~= encBits);
fprintf('📉 BER (öncesi): %.5f\n', ber_pre);

%% 11. LDPC decode
decBits = zeros(numBl*K,1);
for i=1:numBl
    decBits((i-1)*K+1:i*K) = ldpcDecode(llr((i-1)*N+1:i*N), decCfg, 50);
end

%% 12. CRC kontrol
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

%% 13. output_1.txt (tüm alınan ham bitler)
write_bits_to_file('output_1.txt', rxHard);

%% 14. output_2.txt (CRC sonrası alınan ham bitler)
write_bits_to_file('output_2.txt', rxCRC);

%% 15. Konstelasyon diyagramı
figure;
scatter(real(symbols(1:min(1000,end))), imag(symbols(1:min(1000,end))), '.b'); hold on;
scatter(real(rxSym(1:min(1000,end))), imag(rxSym(1:min(1000,end))), '.r');
title(['Konstelasyon: ',modType,' @',num2str(EbNo_dB),' dB']); legend('Tx','Rx');

%% 16. LLR histogramı
figure; histogram(llr,100); title('LLR Dağılımı');

%% 17. İlk 100 bit karşılaştırması (altlı üstlü)
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
