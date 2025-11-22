clc; clear; close all;

%% GENEL PARAMETRELER
fs = 1e6;                  % Örnekleme frekansı (Hz)
T = 1;                     % Sinyal süresi (sn)
t = 0:1/fs:T-1/fs;         % Zaman vektörü

center_freqs = [100e3 250e3 400e3]; % Kullanılabilir merkez frekanslar (Hz)
jammer_freqs = [250e3];     % Birden fazla jammer için [250e3, 400e3] gibi ekleyebilirsin
jammer_amp = 2.5;           % Jammer genliği
noise_level = 0.1;          % Gürültü genliği

%% BAŞLANGIÇ FREKANSI (İlk olarak jammed olan frekansta başlatıyoruz)
current_center = 250e3;
fprintf("Başlangıç merkez frekansı: %.1f kHz\n", current_center/1e3);

%% SİNYAL ÜRETİMİ (Ana sinyal + tüm jammer'lar + gürültü)
signal = sin(2*pi*current_center*t);
jammer = zeros(size(t));
for jf = jammer_freqs
    jammer = jammer + jammer_amp*sin(2*pi*jf*t);
end
noise = noise_level*randn(size(t));
rx = signal + jammer + noise;

%% PSD ANALİZİ
N = 2^nextpow2(length(rx));
f = fs*(0:(N/2))/N;
RX = fft(rx, N);
PSD = abs(RX/N).^2;
PSD = PSD(1:N/2+1);

% Grafik: Spektrumun genel görünümü
figure; plot(f/1e3,10*log10(PSD),'b','LineWidth',1.2);
xlabel('Frekans (kHz)'); ylabel('Güç (dB)');
title('Spektrum (PSD)');
grid on;

%% JAMMER DETECTION (EN GELİŞMİŞ YÖNTEM)
% Her bir merkez frekans için:
% - Local background/median alınır
% - Training/guard cells ile local threshold belirlenir
% - Sadece merkez çevresinde karar verilir

window = 4000;     % Her merkez frekansta +- 2 kHz'lik aralığı inceleriz
training_width = 10000; % Training cell aralığı (daha geniş arka plan ölçümü için, +-5 kHz)
threshold_dB = 10; % Local medyandan 10 dB üstü jammer olarak sayılır

jammed_centers = [];
clean_centers = [];

for i = 1:length(center_freqs)
    % Frekans aralığı indeksleri
    idx_center = find(abs(f - center_freqs(i)) < window);    % Sinyal+jammer aralığı
    idx_training = find(abs(f - center_freqs(i)) > window & abs(f - center_freqs(i)) < training_width);
    
    % Local background (training cell) dB olarak medyanını al (robust, outlier etkisiz)
    local_bg = median(10*log10(PSD(idx_training)));
    % O frekansta max güç (ana sinyal + varsa jammer dâhil)
    local_sig = max(10*log10(PSD(idx_center)));
    
    % Karar: Sinyal gücü local background'ın threshold_dB üstünde mi?
    if (local_sig - local_bg) > threshold_dB
        jammed_centers = [jammed_centers center_freqs(i)];
    else
        clean_centers = [clean_centers center_freqs(i)];
    end
end

if isempty(jammed_centers)
    disp('Hiçbir merkez frekansta jammer yok! İletişime devam.');
else
    fprintf('Jammer tespit edilen merkez frekans(lar): %.1f kHz\n', jammed_centers/1e3);
end

fprintf('Temiz kalan merkez frekans(lar): %.1f kHz\n', clean_centers/1e3);

%% KAÇIŞ (FREKANS ATLAMA) & YENİ PSD
if isempty(clean_centers)
    disp('Tüm frekanslar jammed! Kaçacak yer yok.');
else
    next_center = clean_centers(1);
    fprintf('Frekans atlanıyor! Yeni merkez frekans: %.1f kHz\n', next_center/1e3);

    % Temiz frekansta iletim-alım (ve varsa yine jammer)
    signal2 = sin(2*pi*next_center*t);
    jammer2 = zeros(size(t));
    for jf = jammer_freqs
        jammer2 = jammer2 + jammer_amp*sin(2*pi*jf*t);
    end
    rx_new = signal2 + jammer2 + noise;

    RX2 = fft(rx_new, N);
    PSD2 = abs(RX2/N).^2;
    PSD2 = PSD2(1:N/2+1);

    figure;
    plot(f/1e3,10*log10(PSD2),'g','LineWidth',1.2);
    hold on;
    % Eğer hâlâ başka jammer varsa (ör: birden fazla jam varsa)
    for i = 1:length(jammer_freqs)
        xline(jammer_freqs(i)/1e3, '--r', ['Jam ', num2str(jammer_freqs(i)/1e3),' kHz']);
    end
    xlabel('Frekans (kHz)'), ylabel('Güç (dB)')
    title(['Yeni Frekansta PSD (', num2str(next_center/1e3),' kHz)'])
    legend('PSD (Frekans Atlandı)','Jam. Frekansları')
    grid on

    disp('Sistem yeni frekansta temiz şekilde çalışıyor.');
end
