function adapt_mod_ldpc_crc_jammer_video()
clc; close all;

%% 1. Kullanıcıdan SNR ve sistem parametreleri
EbNo_dB = input('SNR (dB) değerini girin: ');
center_freqs = [100e3 250e3 400e3];
jammer_freq = 250e3;
jammer_amp = 2;
noise_level = 0.15;
fs = 1e6;
T = 1;
t = 0:1/fs:T-1/fs;

%% 2. Video dosyasını oku
fid = fopen("C:\Users\mfb36\Desktop\aa-MediComm\aa_MediComm_Sim_Files\MediComm_Sim_Codes\MATLAB\input_video\input_video.mp4", 'rb');
if fid == -1
    error('video.mp4 bulunamadı!');
end
video_bytes = fread(fid, 'uint8');
fclose(fid);

chunk_size = 100000;
total_bytes = length(video_bytes);
num_chunks = ceil(total_bytes / chunk_size);
fprintf("Toplam parça sayısı: %d\n", num_chunks);

all_rx_bytes = [];

%% LDPC & CRC sabit yapılandırmaları
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
K = encCfg.NumInformationBits;
N = encCfg.BlockLength;
crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');

for part = 1:num_chunks
    fprintf("\n🧩 Parça %d/%d işleniyor...\n", part, num_chunks);

    % Parçayı al
    idx_start = (part-1)*chunk_size + 1;
    idx_end = min(part*chunk_size, total_bytes);
    part_bytes = video_bytes(idx_start:idx_end);

    % Byte → Bit
    video_bits = reshape(de2bi(part_bytes, 8, 'left-msb')', [], 1);

    %% 3. CRC ekle
    txCRC = crcGen(video_bits);

    %% 4. Padding
    pad = mod(K - mod(length(txCRC),K),K);
    if pad > 0
        txCRC = [txCRC; zeros(pad,1)];
    end

    %% 5. LDPC encode
    numBl = length(txCRC)/K;
    encBits = zeros(numBl*N,1);
    for i = 1:numBl
        encBits((i-1)*N+1:i*N) = ldpcEncode(txCRC((i-1)*K+1:i*K), encCfg);
    end

    %% 6. JAMMER TESPİTİ + FREKANS ATLAMA
    fprintf('Frekans taraması başlıyor...\n');
    clean_centers = [];
    for cf = center_freqs
        if abs(cf - jammer_freq) < 1
            jammer = jammer_amp*sin(2*pi*jammer_freq*t);
        else
            jammer = zeros(size(t));
        end
        noise = noise_level*randn(size(t));
        rx_total = jammer + noise;

        Nfft = 2^nextpow2(length(rx_total));
        f = fs*(0:(Nfft/2))/Nfft;
        RX = fft(rx_total, Nfft);
        PSD = abs(RX/Nfft).^2;
        PSD = PSD(1:Nfft/2+1);

        margin = 500;
        idx = find(abs(f - cf) < margin);
        [~, sortidx] = sort(PSD(idx), 'descend');
        top_vals = PSD(idx(sortidx(1:min(5,end))));
        power_dB = 10*log10(mean(top_vals));
        bg_band = PSD(1:round(length(PSD)*0.25));
        bg_dB = 10*log10(mean(bg_band));
        threshold = bg_dB + 15;

        fprintf('%.0f kHz: Güç = %.2f dB, BG = %.2f dB, Eşik = %.2f dB -> ', cf/1e3, power_dB, bg_dB, threshold);
        if power_dB > threshold
            fprintf('JAMMED\n');
        else
            fprintf('TEMİZ\n');
            clean_centers = [clean_centers, cf];
        end
    end

    %% 7. Temiz frekans seçimi
    if isempty(clean_centers)
        error('Kaçacak frekans kalmadı.');
    else
        new_center = clean_centers(1);
        fprintf('\nKullanılacak frekans: %.0f kHz\n', new_center/1e3);
    end

    %% 8. Kanal ortamı
    signal = sin(2*pi*new_center*t);
    if abs(new_center - jammer_freq) < 1
        jammer = jammer_amp*sin(2*pi*jammer_freq*t);
    else
        jammer = zeros(size(t));
    end
    noise = noise_level*randn(size(t));
    rx_total = signal + jammer + noise;

    %% 9. Modülasyon
    modType = select_modulation(EbNo_dB);
    fprintf('Seçilen modülasyon: %s\n', modType);

    switch modType
        case 'BPSK'
            L = 1; M = 2;
            symbols = 1 - 2*double(encBits);
            Eb = 1;
            N0 = Eb / (10^(EbNo_dB/10));
            noise_ch = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
            rxSym = symbols + noise_ch;
            llr = (2/N0)*real(rxSym);
        otherwise
            L = log2(str2double(modType(1:end-3)));
            M = 2^L;
            bits = reshape(encBits, L, []).';
            symbols = qammod(bi2de(bits, 'left-msb'), M, 'UnitAveragePower', true);
            Es = mean(abs(symbols).^2);
            Eb = Es / L;
            N0 = Eb / (10^(EbNo_dB/10));
            noise_ch = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
            rxSym = symbols + noise_ch;
            llr = qamdemod(rxSym, M, 'OutputType', 'approxllr', ...
                'UnitAveragePower', true, 'NoiseVariance', N0);
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
    [rxCRC, err] = crcDet(decBits);
    ber_post = mean(decBits(1:length(txCRC)) ~= txCRC);
    fprintf('📈 BER (sonrası): %.5f\n', ber_post);

    if err
        warning('❌ CRC başarısız!');
        rxCRC = decBits(1:end - pad);
    else
        disp('✅ CRC başarılı.');
        if pad > 0 && length(rxCRC) >= pad
            rxCRC = rxCRC(1:end-pad);
        end
    end

    %% 13. Bit → Byte ve biriktirme
    padbits = mod(8 - mod(length(rxCRC), 8), 8);
    rxCRC = [rxCRC; zeros(padbits,1)];
    video_bytes_rx = bi2de(reshape(rxCRC, 8, [])', 'left-msb');
    all_rx_bytes = [all_rx_bytes; video_bytes_rx];
end

%% 14. Video dosyasına geri yaz
fid = fopen("C:\Users\mfb36\Desktop\aa-MediComm\aa_MediComm_Sim_Files\MediComm_Sim_Codes\MATLAB\output_video\received_video.mp4", 'wb');
fwrite(fid, all_rx_bytes, 'uint8');
fclose(fid);
fprintf('🎥 received_video.mp4 dosyası yazıldı.\n');

%% 15. Konstelasyon diyagramı
figure;
scatter(real(symbols(1:min(1000,end))), imag(symbols(1:min(1000,end))), '.b'); hold on;
scatter(real(rxSym(1:min(1000,end))), imag(rxSym(1:min(1000,end))), '.r');
title(['Konstelasyon: ',modType,' @',num2str(EbNo_dB),' dB']); legend('Tx','Rx');

%% 16. LLR histogramı
figure; histogram(llr,100); title('LLR Dağılımı');

end

function m = select_modulation(EbNo_dB)
    if EbNo_dB < 4, m = 'BPSK';
    elseif EbNo_dB < 8, m = '4QAM';
    elseif EbNo_dB < 15, m = '16QAM';
    else, m = '64QAM';
    end
end
