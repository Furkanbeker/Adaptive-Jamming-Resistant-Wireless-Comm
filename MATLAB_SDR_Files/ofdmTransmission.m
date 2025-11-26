function [legitimacy, txSerial] = ofdmTransmission()
% ofdmTransmission: Pluto TX->RX sabit-pilotlu OFDM çerçevesi gönder/al,
% senkronizasyon + CFO düzeltme yap, kanal kestirimi çıkar ve görsele kaydet.
%
% Çıktılar:
%  - legitimacy : 'legitimate' | 'illegitimate' | 'unknown' (basit kural)
%  - txSerial   : TX Pluto seri numarası
%
% Notlar:
%  - Bu tek dosyadır; m-sequence için local_mseq alt fonksiyonu kullanılır.
%  - Communications Toolbox Support Package for ADALM-PLUTO gereklidir.

clc;
clearvars -except legitimacy txSerial;

%% ---------- Pluto bul ve eşle ----------
radios = findPlutoRadio; % Support Package fonksiyonu
if numel(radios) < 2
    error('En az 2 adet ADALM-Pluto bağlı olmalı. Bulunan: %d', numel(radios));
end
txSerial = radios(1).SerialNum;
rxSerial = radios(2).SerialNum;

hw.Tx_ID = ['sn:' num2str(txSerial)];
hw.Rx_ID = ['sn:' num2str(rxSerial)];

%% ---------- OFDM parametreleri ----------
FFTsize              = 64;
NumActiveSubcarriers = 48;                 % merkezde 48 aktif taşıyıcı
CPsize               = 8;
NumSymbols           = 16;                 % çerçeve başına 16 OFDM sembol
df                   = 15e3;               %#ok<NASGU> alt taşıyıcı aralığı (bilgi amaçlı)

% Aktif taşıyıcı indeksleri (1-bazlı; DC merkeze gelecek şekilde)
active_carriers = (FFTsize - NumActiveSubcarriers)/2 + (1:NumActiveSubcarriers);
assert(mod(FFTsize,2)==0, 'FFTsize çift olmalı.');

%% ---------- Modülasyon (pilot-only çerçeve) ----------
Mmod   = 16;                 % 16-QAM
Mb     = log2(Mmod);
nAct   = NumActiveSubcarriers;

% Sabit pilot sembolleri (aktif taşıyıcı sayısı kadar) – tüm sembollerde aynı
bits   = randi([0 1], nAct*Mb, 1);
symIdx = bi2de(reshape(bits, Mb, []).','left-msb');
pilotF = qammod(symIdx, Mmod, 'UnitAveragePower', true);   % nAct x 1

% Frekans alanı OFDM matrisi: (FFTsize x NumSymbols)
OFDM_F = zeros(FFTsize, NumSymbols);
OFDM_F(active_carriers, :) = repmat(pilotF, 1, NumSymbols);

% Zaman alanı (IFFT kolonsal)
ofdm_t = ifft(OFDM_F, FFTsize, 1);

% CP ekle ( (FFTsize+CPsize) x NumSymbols )
ofdm_cp = [ofdm_t(end-CPsize+1:end, :); ofdm_t];

% Seri hale getir
tx_ofdm = ofdm_cp(:);
tx_ofdm = tx_ofdm ./ max(abs(tx_ofdm)+eps);

%% ---------- Preamble (M-sequence, 2^8-1 = 255, iki kez) ----------
preamble = [local_mseq(8, 1); local_mseq(8, 1)];  % ±1 kolon vektör

%% ---------- RRC yükseltme filtrelemesi ----------
rrc.alpha = 0.5;
rrc.span  = 12;   % sembol
rrc.sps   = 8;    % upsample oranı

rrc.h = rcosdesign(rrc.alpha, rrc.span, rrc.sps, 'normal');

% Başına durgunluk (DC transiyent azaltma)
padZeros = 500;
baseband = [zeros(padZeros,1); preamble; tx_ofdm];

% Upsample + pulse shaping
up = upsample(baseband, rrc.sps);
tx_frame = conv(up, rrc.h, 'full');

%% ---------- Pluto RF ayarları ----------
hw.Tx_G = -10;          % dB
hw.Rx_G =  25;          % dB
hw.Fc   = 1.00e9;       % Hz
hw.Fs   = 20e6;         % Hz (örnekleme)
rxN     = 3 * length(tx_frame);   % tek tampon tahsisi

Tx = sdrtx('Pluto', ...
    'Gain', hw.Tx_G, ...
    'CenterFrequency', hw.Fc, ...
    'BasebandSampleRate', hw.Fs, ...
    'RadioID', hw.Tx_ID);

Rx = sdrrx('Pluto', ...
    'SamplesPerFrame', rxN, ...
    'OutputDataType','double', ...
    'CenterFrequency', hw.Fc, ...
    'BasebandSampleRate', hw.Fs, ...
    'GainSource','manual', ...
    'Gain', hw.Rx_G, ...
    'RadioID', hw.Rx_ID);

% Sürekli döndür
Tx.transmitRepeat(tx_frame);

% Tek tampon oku (gerekirse arttırılabilir)
rx_raw = Rx();
rx_raw = rx_raw ./ max(abs(rx_raw)+eps);

%% ---------- Kaba zaman tespiti (güç kenarı) ----------
W = 600;                                 % pencere
pow = abs(rx_raw).^2;
% kayan toplam + fark (yükselen kenar arama)
movsum = conv(pow, ones(W,1), 'valid');
dd = [0; diff(movsum)];
[~, idxTop] = maxk(dd, 3);
idxTop = sort(idxTop);
start_coarse = idxTop(min(2,numel(idxTop)));  % 2. en büyük kenarı al
start_coarse = max(start_coarse, 1);

% Preamble aralığı (downsample öncesi)
preLenSamp = numel(preamble) * rrc.sps;
stop_coarse = min(start_coarse + preLenSamp - 1, numel(rx_raw));

seg = rx_raw(start_coarse:stop_coarse);

%% ---------- CFO kestirimi ve düzeltme ----------
% BPSK-benzeri için kare alma CFO tepeyi 2x gösterir → /2
seg2 = seg.^2;
Nfft_hi = 100 * length(seg2);          % ince frekans ızgarası
S = abs(fftshift(fft(seg2, Nfft_hi)));
df = hw.Fs / Nfft_hi;
freq = (-hw.Fs/2 : df : (hw.Fs/2 - df)).';
[~, kmax] = max(S);
f_off = freq(kmax) / 2;                % /2 kritik

t = (0:numel(rx_raw)-1).' / hw.Fs;
rx_cfo = rx_raw .* exp(-1j*2*pi*f_off*t);

%% ---------- Matched filter + optimum örnekleme fazı ----------
rx_mf = conv(rx_cfo, rrc.h, 'full');

% Grup gecikmesini telafi et (RRC toplam gecikme: span*sps/2 örnek)
gd = (rrc.span * rrc.sps)/2;
rx_mf = rx_mf( gd+1 : gd+numel(rx_cfo) );  % uzunluğu rx_cfo ile hizala

% Faz tarama (en yüksek güç hangi fazda?)
p = zeros(rrc.sps,1);
for k = 1:rrc.sps
    tmp = rx_mf(k:rrc.sps:end);
    p(k) = mean(abs(tmp).^2);
end
[~, bestPhase] = max(p);
rx_ds = rx_mf(bestPhase:rrc.sps:end);     % downsample edilmiş

%% ---------- İnce zaman (preamble korelasyonu) ----------
% dwn ile preamble korelasyonu: lag -> başlangıç
[c, lags] = xcorr(rx_ds, preamble);
[~, imx]  = max(abs(c));
lag       = lags(imx);
idx0      = lag + 1;                      % rx_ds indeksinde preamble başlangıcı
idx0      = max(1, idx0);

% Faz düzeltmesi (global)
Lpre = numel(preamble);
if idx0 + Lpre - 1 > numel(rx_ds)
    warning('Preambleyı tam çıkaracak kadar örnek yok.');
    legitimacy = 'unknown';
    release(Rx); release(Tx);
    return;
end
rx_pre = rx_ds(idx0 : idx0+Lpre-1);
phi = angle(sum(conj(preamble) .* rx_pre));
rx_ds = rx_ds .* exp(-1j*phi);

% Preamble sonrası OFDM blok
ofdm_ser_len = numel(tx_ofdm);  % beklenen uzunluk
i1 = idx0 + Lpre;
i2 = i1 + ofdm_ser_len - 1;
if i2 > numel(rx_ds)
    warning('OFDM verisinin tamamı alınamadı (kısa çerçeve).');
    legitimacy = 'unknown';
    release(Rx); release(Tx);
    return;
end
rx_ofdm_ser = rx_ds(i1:i2);

% Matrise dök: (FFTsize+CPsize) x NumSymbols
rx_ofdm_mtx = reshape(rx_ofdm_ser, FFTsize+CPsize, NumSymbols);

% CP at → FFT (kolon bazlı)
rx_fft = fft( rx_ofdm_mtx(CPsize+1:end, :), FFTsize, 1 );

%% ---------- Kanal kestirimi (aktif taşıyıcılar) ----------
H_est = rx_fft(active_carriers, :) ./ (OFDM_F(active_carriers, :) + eps);

% Görsel için normalize et (0..1), merkez 48’e denk gelen 48xNumSymbols zaten
imgM = abs(H_est);
imgM = imgM ./ max(imgM(:) + eps);

% Klasör ve dosya
outDir = 'data_ofdm';
if ~exist(outDir, 'dir'), mkdir(outDir); end
iter = 1;
fname = fullfile(outDir, sprintf('1_iter_%d.jpg', iter));  % başına '1_' kondu
imwrite(imgM, fname, 'jpg');

% Basit sınıflandırma kuralı (dosya adına göre değil, sinyal niteliğine göre)
avgGain = mean(imgM(:));
if avgGain > 0.35
    legitimacy = 'legitimate';
elseif avgGain < 0.15
    legitimacy = 'illegitimate';
else
    legitimacy = 'unknown';
end

% Kaynakları bırak
release(Rx);
release(Tx);

end % main function


%% ====== Yardımcı: M-sequence (LFSR) üretici (±1, uzunluk 2^m-1) ======
function ms = local_mseq(m, shift)
% local_mseq(m, shift): m-dereceli maksimum uzunluk sekansı (±1).
% Varsayılan polinom: x^8 + x^4 + x^3 + x^2 + 1 (m=8 için)
% Giriş:
%   m     : LFSR derecesi (ör. 8)
%   shift : opsiyonel, döngüsel kaydırma (ör. 1)
% Çıkış:
%   ms    : (2^m-1)x1, değerler ±1 (double)

if nargin < 2, shift = 0; end
N = 2^m - 1;

% m=8 için yaygın primitive taps:
% Polinom: x^8 + x^4 + x^3 + x^2 + 1 → tap indexleri [8 4 3 2]
% Genel durumda tabloya göre taps seçmek gerekir; burada m=8 hedeflendi.
if m == 8
    taps = [8 4 3 2];
else
    error('local_mseq şu an m=8 için tanımlı. İhtiyaç varsa taps tablosu eklenmeli.');
end

reg = ones(1, m);     % başlangıç koşulu (hepsi 1)
out = zeros(N,1);

for i = 1:N
    out(i) = reg(end);
    fb = mod(sum(reg(taps)), 2);     % feedback biti
    reg = [fb reg(1:end-1)];
end

% {0,1} → {+1, -1} (BPSK haritalama)
ms = 2*out - 1;

% Döngüsel kaydırma (ileri yönde)
if shift ~= 0
    ms = circshift(ms, shift);
end
end
