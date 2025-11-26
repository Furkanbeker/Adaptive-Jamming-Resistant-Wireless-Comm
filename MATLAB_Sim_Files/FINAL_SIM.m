function ldpc_jammer_gui_allinone()
% MediComm 5-Band Live (RF=2.4 GHz, IF simülasyon) — Rev5
% - 5 bant RF: 2402/2406/2410/2414/2418 MHz (UI'da GHz)
% - IF simülasyon: 100/200/300/400/500 kHz, fs=1.2 MHz (500 kHz Nyquist güvenli)
% - Jammer 5 buton, Jammer otomatik gezme (opsiyon)
% - SNR slider + "SNR otomatik rastgele" (her parçada yeni SNR)
% - Sol altta: Dosya Boyutu / İletilen / Kalan / Hız / ETA (canlı)
% - İlerleme: linear gauge
% - Spektrum: çizgiler ince & yarı saydam; bant çizgileri ince gri
% - YENİ: DURDUR -> DEVAM ET (gri) ve kaldığı yerden devam et

%% -------------------- SABİTLER (ŞARTNAME) --------------------
BANDS_RF_HZ = [2.402e9, 2.406e9, 2.410e9, 2.414e9, 2.418e9];   % RF merkezler
BANDS_IF_HZ = [100e3, 200e3, 300e3, 400e3, 500e3];            % IF merkezler (sim)

%% -------------------- STATE --------------------
state.running       = false;     % döngü içinde mi
state.paused        = false;     % duraklatıldı mı
state.stopRequested = false;     % döngüde çıkış işareti
state.simInit       = false;     % simülasyon ilk kez hazırlandı mı
state.videoPath     = "";
state.outputFolder  = "";
state.fs            = 1.2e6;     % IF örnekleme (Nyquist 600 kHz)
state.T             = 1;         % 1 sn pencere (sadece spektrum/algı için)
state.bandsRF       = BANDS_RF_HZ(:).';
state.bandsIF       = BANDS_IF_HZ(:).';
state.jammer_amp    = 2.0;
state.EbNo_dB       = 6;
state.chunk_size    = 100000;
state.maxItersLDPC  = 50;
state.progress      = 0;
state.all_rx_bytes  = uint8([]);
state.jamAutoWander = false;
state.jamIdx        = 3;         % jammer bandı
state.cfIdx         = 3;         % çalışma bandı
state.snrAuto       = false;     % SNR otomatik

% çalışma oturumu (resume için gerekli)
state.sim = struct('video_bytes',[], 'total_bytes',0, 'num_chunks',0, ...
                   'outFile',"", 'part_idx',0, 'elapsedActive',0, ...
                   'berPreSeries',[], 'berPostSeries',[]);

% LDPC/CRC
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
state.K = encCfg.NumInformationBits;
state.N = encCfg.BlockLength;
state.crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
state.crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');

%% -------------------- UI --------------------
fig = uifigure('Name','MediComm 5-Band Live','Position',[60 60 1420 860]);
gl  = uigridlayout(fig,[3,3], 'RowHeight',{'fit','1x','fit'}, ...
    'ColumnWidth',{400,'1x',560}, 'Padding',8,'RowSpacing',8,'ColumnSpacing',10);

% -------- Sol panel --------
pLeft = uipanel(gl,'Title','Kontroller','FontWeight','bold'); pLeft.Layout.Row=[1 3]; pLeft.Layout.Column=1;
gL = uigridlayout(pLeft,[23,2], 'RowHeight',repmat({'fit'},1,23), ...
    'ColumnWidth',{140,240}, 'RowSpacing',6,'ColumnSpacing',10,'Padding',8);

% Video/Çıkış
uibutton(gL,'Text','Video Seç (mp4/bin)','ButtonPushedFcn',@(btn,~) onChooseVideo());
uilabel(gL,'Text',' ');

uilabel(gL,'Text','Çıkış klasörü:','HorizontalAlignment','right');
uibutton(gL,'Text','Seç','ButtonPushedFcn',@(btn,~) onChooseOutput());
lblOut = uilabel(gL,'Text','(varsayılan: girişle aynı klasör)','WordWrap','on'); lblOut.Layout.Column=[1 2];

% SNR
uilabel(gL,'Text','SNR (dB):','HorizontalAlignment','right');
sSNR = uislider(gL,'Limits',[-2 20],'Value',state.EbNo_dB, 'MajorTicks',-2:2:20,'MinorTicks',[], ...
    'ValueChangedFcn',@(s,~) setSNR(s.Value));
sSNR.Layout.Column=[1 2];

uilabel(gL,'Text','SNR otomatik:','HorizontalAlignment','right');
cbSNRAuto = uicheckbox(gL,'Text','Rastgele değişsin','Value',state.snrAuto,'ValueChangedFcn',@(c,~) setSNRAuto(c.Value));

% Jammer
uilabel(gL,'Text','Jammer amp:','HorizontalAlignment','right');
sJamA = uislider(gL,'Limits',[0 5],'Value',state.jammer_amp,'MajorTicks',0:1:5,'MinorTicks',[],...
    'ValueChangedFcn',@(s,~) setJamAmp(s.Value)); sJamA.Layout.Column=[1 2];

uilabel(gL,'Text','Jammer Bandı:','HorizontalAlignment','right');
jamBtnPanel = uipanel(gL,'BorderType','none'); jamBtnPanel.Layout.Column=2;
gJB = uigridlayout(jamBtnPanel,[1,5],'ColumnWidth',{'1x','1x','1x','1x','1x'},'RowHeight',{'fit'},'ColumnSpacing',6);
btnJ = gobjects(1,5);
for k=1:5
    btnJ(k) = uibutton(gJB,'Text',sprintf('B%d (%.3f)',k,state.bandsRF(k)/1e9), 'ButtonPushedFcn',@(b,~) selectJammer(k));
end

uilabel(gL,'Text','Jammer otomatik:','HorizontalAlignment','right');
cbWander = uicheckbox(gL,'Text','Band gezer','Value',state.jamAutoWander,'ValueChangedFcn',@(c,~) setWander(c.Value));

% Chunk
uilabel(gL,'Text','Chunk size (B):','HorizontalAlignment','right');
eChunk = uieditfield(gL,'numeric','Limits',[1 Inf],'Value',state.chunk_size,'RoundFractionalValues','on','ValueDisplayFormat','%.0f',...
    'ValueChangedFcn',@(e,~) setChunkSize(e.Value));

% Başlat / Durdur-Devam
uilabel(gL,'Text',' ','HorizontalAlignment','right'); % spacer
pStart = uipanel(gL,'BorderType','none'); pStart.Layout.Column=2;
gStart = uigridlayout(pStart,[1,2],'ColumnWidth',{120,120},'ColumnSpacing',8);
btnStart = uibutton(gStart,'Text','BAŞLAT','BackgroundColor',[0.21 0.6 0.24], 'FontWeight','bold','ButtonPushedFcn',@(b,~) onStartNew());
btnStop  = uibutton(gStart,'Text','DURDUR','BackgroundColor',[0.75 0.2 0.2], 'Enable','off','ButtonPushedFcn',@(b,~) onStopOrResume());

% Seçili video / durum
uilabel(gL,'Text','Seçili video:','HorizontalAlignment','right'); lblVideo = uilabel(gL,'Text','yok','WordWrap','on');
uilabel(gL,'Text','Durum:','HorizontalAlignment','right');      lblStats = uilabel(gL,'Text','hazır','WordWrap','on');

% İlerleme
uilabel(gL,'Text','İlerleme:','HorizontalAlignment','right'); gauge = uigauge(gL,'linear','Limits',[0 1],'Value',0); gauge.Layout.Column=2;
uilabel(gL,'Text',' ','HorizontalAlignment','right'); lblProg = uilabel(gL,'Text','0 %');

% Aktarım Durumu
pXfer = uipanel(gL,'Title','Aktarım Durumu','FontWeight','bold'); pXfer.Layout.Column=[1 2];
gx = uigridlayout(pXfer,[5,2],'RowHeight',repmat({'fit'},1,5),'ColumnWidth',{140,220},'RowSpacing',4,'Padding',6,'ColumnSpacing',10);
uilabel(gx,'Text','Dosya Boyutu:','HorizontalAlignment','right');  lblSize   = uilabel(gx,'Text','-');
uilabel(gx,'Text','İletilen:','HorizontalAlignment','right');      lblSent   = uilabel(gx,'Text','-');
uilabel(gx,'Text','Kalan:','HorizontalAlignment','right');         lblRemain = uilabel(gx,'Text','-');
uilabel(gx,'Text','Gönderme Hızı:','HorizontalAlignment','right'); lblRate   = uilabel(gx,'Text','-');
uilabel(gx,'Text','Kalan Süre:','HorizontalAlignment','right');    lblETA    = uilabel(gx,'Text','-');

% -------- Orta: Grafikler --------
pMid = uipanel(gl,'Title','Grafikler','FontWeight','bold'); pMid.Layout.Row=[1 2]; pMid.Layout.Column=2;
gM = uigridlayout(pMid,[3,1],'RowHeight',{'1x','1x','1x'},'Padding',6,'RowSpacing',8);
axSpec  = uiaxes(gM); title(axSpec,'Spektrum (IF, tıkla: en yakın banda snap)'); xlabel(axSpec,'Frekans (IF, Hz)'); ylabel(axSpec,'PSD');
axSpec.HitTest='on'; axSpec.PickableParts='all'; axSpec.ButtonDownFcn=@(~,~) onSpectrumClick();
axConst = uiaxes(gM); title(axConst,'Konstelasyon'); xlabel(axConst,'I'); ylabel(axConst,'Q');
axLLR   = uiaxes(gM); title(axLLR,'LLR Histogramı'); xlabel(axLLR,'LLR'); ylabel(axLLR,'Adet');

% -------- Sağ: Metrikler --------
pRight = uipanel(gl,'Title','Metrikler & Kayıt','FontWeight','bold'); pRight.Layout.Row=[1 2]; pRight.Layout.Column=3;
gR = uigridlayout(pRight,[6,1],'RowHeight',{'fit','fit','1x','fit','1x','1x'},'Padding',8,'RowSpacing',10);

gTop = uigridlayout(gR,[3,2],'RowHeight',{'fit','fit','fit'},'ColumnWidth',{160,160},'RowSpacing',6,'ColumnSpacing',10);
uilabel(gTop,'Text','SNR (dB):','HorizontalAlignment','right');         lblSNR = uilabel(gTop,'Text',num2str(state.EbNo_dB,'%0.2f'));
uilabel(gTop,'Text','Çalışma CF (GHz):','HorizontalAlignment','right'); lblCF  = uilabel(gTop,'Text',num2str(state.bandsRF(state.cfIdx)/1e9,'%.3f'));
uilabel(gTop,'Text','Jammer f (GHz):','HorizontalAlignment','right');   lblJam = uilabel(gTop,'Text',num2str(state.bandsRF(state.jamIdx)/1e9,'%.3f'));

gBER = uigridlayout(gR,[2,2],'RowHeight',{'fit','fit'},'ColumnWidth',{160,160});
uilabel(gBER,'Text','BER (öncesi):','HorizontalAlignment','right');  lblBERpre  = uilabel(gBER,'Text','-');
uilabel(gBER,'Text','BER (sonrası):','HorizontalAlignment','right'); lblBERpost = uilabel(gBER,'Text','-');

axBER = uiaxes(gR); title(axBER,'BER Zaman Serisi'); xlabel(axBER,'Parça #'); ylabel(axBER,'BER');
hold(axBER,'on'); pltPre = plot(axBER,nan,nan,'-'); pltPost = plot(axBER,nan,nan,'-'); legend(axBER,{'Pre-LDPC','Post-LDPC'},'Location','northeast');

uilabel(gR,'Text','5 Bant Taraması');
tbl = uitable(gR,'ColumnName',{'Band','RF_GHz','Power_dB','BG_dB','Thresh','Durum'},'Data',cell(0,6));
ta  = uitextarea(gR,'Editable','off','FontName','Consolas');

lblHelp = uilabel(gl,'Text','İpucu: Jammer butonları RF (GHz). Spektrumda tıklayınca IF’te en yakın banda snap eder; metrikler RF gösterir.');
lblHelp.Layout.Row=3; lblHelp.Layout.Column=[2 3];

%% -------------------- CALLBACKS --------------------
    function onChooseVideo()
        [f,p] = uigetfile({'*.mp4;*.bin;*.dat','Video/Raw (*.mp4,*.bin,*.dat)';'*.*','Hepsi'});
        if isequal(f,0), return; end
        state.videoPath = fullfile(p,f);
        lblVideo.Text = " " + f;
        logmsg("Video seçildi: " + state.videoPath);
    end

    function onChooseOutput()
        p = uigetdir();
        if isequal(p,0), return; end
        state.outputFolder = p;
        lblOut.Text = p;
        logmsg("Çıkış klasörü: " + p);
    end

    function setSNR(v), state.EbNo_dB = v; lblSNR.Text = sprintf('%.2f',v); end
    function setSNRAuto(tf), state.snrAuto = tf; sSNR.Enable = iff(tf,'off','on'); end
    function setJamAmp(v), state.jammer_amp = v; end
    function setChunkSize(v), v=max(1000,round(v)); state.chunk_size=v; eChunk.Value=v; end
    function setWander(tf)
        state.jamAutoWander = tf;
        for k=1:5, btnJ(k).Enable = iff(tf,'off','on'); end
    end

    function selectJammer(k)
        state.jamIdx = k; lblJam.Text = sprintf('%.3f', state.bandsRF(k)/1e9);
        logmsg(sprintf("Jammer band: B%d (%.3f GHz)",k,state.bandsRF(k)/1e9));
    end

    function onSpectrumClick()
        cp = axSpec.CurrentPoint; x = cp(1,1); % IF Hz
        [~,idx] = min(abs(state.bandsIF - x));
        state.cfIdx = idx; lblCF.Text = sprintf('%.3f',state.bandsRF(idx)/1e9);
        logmsg(sprintf("CF grafikten: B%d (%.3f GHz)",idx,state.bandsRF(idx)/1e9));
    end

    % BAŞLAT: daima yeni oturum başlatır (sıfırdan)
    function onStartNew()
        if state.videoPath == "", uialert(fig,'Lütfen bir video dosyası seçin.','Uyarı'); return; end
        if state.running, return; end

        % Oturum init
        initSimulation();  % state.sim doldurulur
        state.paused = false; state.stopRequested = false; state.running = true;
        btnStart.Enable = 'off';
        btnStop.Enable  = 'on';  btnStop.Text = 'DURDUR'; btnStop.BackgroundColor = [0.75 0.2 0.2];
        lblStats.Text   = 'çalışıyor...';

        % döngü: 1. parçadan
        loopSimulation(state.sim.part_idx + 1);

        % döngü bitti: pause mı finiş mi?
        state.running = false;
        if state.paused
            lblStats.Text = 'durakladı';
            % duraklamada BAŞLAT'ı yeni başlatma için açık bırakmak istersen:
            btnStart.Enable = 'on'; % istenirse 'off' yapıp sadece "DEVAM ET" kullanılabilir
        else
            lblStats.Text = 'bitti';
            btnStart.Enable = 'on'; btnStop.Enable = 'off'; btnStop.Text = 'DURDUR';
        end
    end

    % DURDUR <-> DEVAM ET (toggle)
    function onStopOrResume()
        if state.running
            % DURDUR talebi
            state.stopRequested = true;
            state.paused = true;
            btnStop.Text = 'DEVAM ET'; btnStop.BackgroundColor = [0.6 0.6 0.6];
            lblStats.Text = 'duraklatılıyor...';
        else
            % DEVAM ET
            if ~state.simInit || ~state.paused, return; end
            state.stopRequested = false; state.paused = false; state.running = true;
            btnStop.Text = 'DURDUR'; btnStop.BackgroundColor = [0.75 0.2 0.2];
            btnStart.Enable = 'off'; lblStats.Text = 'çalışıyor...';

            loopSimulation(state.sim.part_idx + 1);

            state.running = false;
            if state.paused
                lblStats.Text = 'durakladı';
            else
                lblStats.Text = 'bitti';
                btnStart.Enable='on'; btnStop.Enable='off'; btnStop.Text='DURDUR';
            end
        end
    end

%% -------------------- ÇEKİRDEK --------------------
    function initSimulation()
        % dosyayı oku ve sim oturumunu hazırla
        [video_bytes,total_bytes] = readBytes(state.videoPath);
        state.all_rx_bytes = uint8([]);
        state.sim.video_bytes = video_bytes;
        state.sim.total_bytes = total_bytes;
        state.sim.num_chunks  = ceil(double(total_bytes)/double(state.chunk_size));
        state.sim.outFile     = resolveOutFile();
        state.sim.part_idx    = 0;          % henüz işlenmemiş
        state.sim.elapsedActive = 0;        % sadece aktif süre toplar
        state.sim.berPreSeries  = nan(1,state.sim.num_chunks);
        state.sim.berPostSeries = nan(1,state.sim.num_chunks);
        state.simInit = true;

        % aktarım kutusu reset
        lblSize.Text   = fmtBytes(total_bytes);
        lblSent.Text   = fmtBytes(0);
        lblRemain.Text = fmtBytes(total_bytes);
        lblRate.Text   = '-'; lblETA.Text = '-';
        gauge.Value    = 0; lblProg.Text = '0 %';

        % sağ BER grafiği reset
        set(pltPre,'XData',nan,'YData',nan); set(pltPost,'XData',nan,'YData',nan);
        tbl.Data = cell(0,6);
    end

    function loopSimulation(start_part)
        % start_part ... num_chunks arasında çalış; duraklama/bitirme durumunda çıkar
        for part = start_part : state.sim.num_chunks
            if state.stopRequested, break; end
            stepTic = tic;

            % --- SNR otomatik rastgele (parça başında) ---
            if state.snrAuto
                snrMin = sSNR.Limits(1); snrMax = sSNR.Limits(2);
                state.EbNo_dB = round(snrMin + (snrMax-snrMin)*rand(), 2);
                sSNR.Value = state.EbNo_dB; lblSNR.Text = sprintf('%.2f', state.EbNo_dB);
            end
            EbNo_dB     = state.EbNo_dB;
            jammer_amp  = state.jammer_amp;

            % Jammer bandı (IF)
            if state.jamAutoWander
                state.jamIdx = randi([1 5]); lblJam.Text = sprintf('%.3f', state.bandsRF(state.jamIdx)/1e9);
            end
            jammer_freq_if = state.bandsIF(state.jamIdx);

            % Çalışma bandı (IF)
            cf_if = state.bandsIF(state.cfIdx);

            % Parça indeksleri
            idx_start = (part-1)*state.chunk_size + 1;
            idx_end   = min(part*state.chunk_size, state.sim.total_bytes);
            part_bytes = state.sim.video_bytes(idx_start:idx_end);

            % --- ENCODE + KANAL + DECODE (kısa) ---
            video_bits = reshape(de2bi(part_bytes,8,'left-msb')',[],1);
            txCRC = state.crcGen(logical(video_bits));
            pad = mod(state.K - mod(length(txCRC),state.K), state.K);
            if pad>0, txCRC = [txCRC; false(pad,1)]; end

            numBl  = length(txCRC)/state.K;
            encBits = false(numBl*state.N,1);
            for i=1:numBl
                encBits((i-1)*state.N+1:i*state.N) = ldpcEncode(txCRC((i-1)*state.K+1:i*state.K),encCfg);
            end

            % 5 bant taraması ve kaçış (IF)
            [clean_if, tableRows] = scanJammer5(jammer_freq_if, jammer_amp, state.bandsIF, state.bandsRF, state.fs, state.T, EbNo_dB);
            tbl.Data = tableRows;
            if ~isempty(clean_if)
                [~,idx] = min(abs(state.bandsIF - clean_if(1)));
                state.cfIdx = idx; cf_if = state.bandsIF(idx);
                lblCF.Text = sprintf('%.3f', state.bandsRF(idx)/1e9);
            end

            % Görsel kanal (IF) -> Spektrum
            t = 0:1/state.fs:state.T-1/state.fs;
            signal = sin(2*pi*cf_if*t);
            jammer = jammer_amp*sin(2*pi*jammer_freq_if*t);
            N0_vis = 1/(10^(EbNo_dB/10)); noise_vis = sqrt(N0_vis/2)*randn(size(t));
            rx_total = signal + jammer + noise_vis;
            plotPSD5(axSpec, rx_total, state.fs, state.bandsIF, state.bandsRF, cf_if, jammer_freq_if);

            % Mod/Demod
            [llr, rxSym, symbols, modType] = modChanDemod(encBits, EbNo_dB); %#ok<ASGLU>
            plotConst(axConst, symbols, rxSym, modType, EbNo_dB);
            plotLLR(axLLR, llr);

            rxHard = llr(:) < 0;
            ber_pre = mean(rxHard ~= encBits);   lblBERpre.Text  = sprintf('%.6f', ber_pre);

            decBits = false(numBl*state.K,1);
            for i=1:numBl
                decBits((i-1)*state.K+1:i*state.K) = ldpcDecode(llr((i-1)*state.N+1:i*state.N), decCfg, state.maxItersLDPC);
            end
            [rxCRC, err] = state.crcDet(decBits); %#ok<NASGU>
            ber_post = mean(decBits(1:length(txCRC)) ~= txCRC); lblBERpost.Text = sprintf('%.6f', ber_post);

            % Çıkış biriktirme (simülasyonda: orijinal parça)
            state.all_rx_bytes = [state.all_rx_bytes; uint8(part_bytes)]; %#ok<AGROW>

            % ----- Aktarım metrikleri -----
            state.sim.part_idx = part;                               % son işlenen
            state.sim.elapsedActive = state.sim.elapsedActive + toc(stepTic);   % sadece aktif süre
            sent_bytes   = double(idx_end);
            remain_bytes = double(state.sim.total_bytes - idx_end);
            rate_Bps = sent_bytes / max(state.sim.elapsedActive,eps);
            rate_bps = rate_Bps * 8;
            eta_s    = remain_bytes / max(rate_Bps,eps);

            lblSize.Text   = fmtBytes(state.sim.total_bytes);
            lblSent.Text   = fmtBytes(sent_bytes);
            lblRemain.Text = fmtBytes(remain_bytes);
            lblRate.Text   = sprintf('%.2f Mbit/s (%.2f MB/s)', rate_bps/1e6, rate_Bps/1e6);
            lblETA.Text    = fmtTime(eta_s);

            % BER serileri & ilerleme
            state.sim.berPreSeries(part)  = ber_pre;
            state.sim.berPostSeries(part) = ber_post;
            set(pltPre,'XData',1:part,'YData',state.sim.berPreSeries(1:part));
            set(pltPost,'XData',1:part,'YData',state.sim.berPostSeries(1:part));
            state.progress = double(idx_end)/double(state.sim.total_bytes);
            updateProgress(state.progress);

            lblStats.Text = sprintf('Parça %d/%d | Mod: %s | CF=%.3f GHz | Jam=%.3f GHz', ...
                                    part,state.sim.num_chunks,modType, ...
                                    state.bandsRF(state.cfIdx)/1e9, state.bandsRF(state.jamIdx)/1e9);
            drawnow limitrate;

            if state.stopRequested, break; end
        end

        % eğer tüm parçalar bittiyse dosyayı yaz
        if state.sim.part_idx >= state.sim.num_chunks && ~isempty(state.all_rx_bytes)
            try
                fid = fopen(state.sim.outFile,'wb'); fwrite(fid,state.all_rx_bytes,'uint8'); fclose(fid);
                logmsg("Çıktı yazıldı: " + state.sim.outFile);
            catch ME
                logmsg("Çıkış yazım hatası: " + ME.message);
            end
        end
    end

%% -------------------- HELPERLAR --------------------
    function outFile = resolveOutFile()
        if state.outputFolder == ""
            [p,~,~] = fileparts(state.videoPath); outFolder = p;
        else
            outFolder = state.outputFolder;
        end
        outFile = fullfile(outFolder,'received_video.mp4');
    end

    function [video_bytes,total_bytes] = readBytes(path)
        try
            fid = fopen(path,'rb'); if fid==-1, error('Dosya açılamadı'); end
            video_bytes = fread(fid,'uint8'); fclose(fid);
            total_bytes = length(video_bytes);
            logmsg(sprintf("Video boyutu: %s", fmtBytes(total_bytes)));
        catch ME %#ok<NASGU>
            logmsg("Dosya okuma hatası. Sentetik test verisi (1 MB)...");
            video_bytes = uint8(randi([0 255],1e6,1)); total_bytes = numel(video_bytes);
        end
    end

    function [clean_centers_if, tableRows] = scanJammer5(jam_f_if, jam_amp, bandsIF, bandsRF, fs, T, EbNo_dB)
        t = 0:1/fs:T-1/fs;
        N0 = 1/(10^(EbNo_dB/10));
        noise = sqrt(N0/2)*randn(size(t));

        clean_centers_if = [];
        tableRows = cell(0,6);
        for k=1:5
            cf_if = bandsIF(k);
            jammer = jam_amp*sin(2*pi*jam_f_if*t);
            rx_total = jammer + noise;

            Nfft = 2^nextpow2(length(rx_total));
            f = fs*(0:(Nfft/2))/Nfft;
            PSD = abs(fft(rx_total,Nfft)/Nfft).^2; PSD = PSD(1:Nfft/2+1);

            margin = 1000; % Hz (1 kHz pencere)
            idxB = find(abs(f - cf_if) < margin);
            if isempty(idxB)
                pwr_dB = -Inf;
            else
                [~,si] = sort(PSD(idxB),'descend');
                top_vals = PSD(idxB(si(1:min(5,numel(idxB)))));
                pwr_dB = 10*log10(mean(top_vals)+eps);
            end
            bg_band = PSD(1:round(numel(PSD)*0.25));
            bg_dB   = 10*log10(mean(bg_band)+eps);
            thr_dB  = bg_dB + 15;

            if pwr_dB > thr_dB
                stat = "JAMMED";
            else
                stat = "TEMIZ";
                clean_centers_if(end+1) = cf_if; %#ok<AGROW>
            end
            tableRows(end+1,:) = {k, bandsRF(k)/1e9, pwr_dB, bg_dB, thr_dB, char(stat)}; %#ok<AGROW>
        end
    end

    function plotPSD5(ax, rx_total, fs, bandsIF, bandsRF, cf_if, jf_if)
        Nfft = 2^nextpow2(length(rx_total));
        f = fs*(0:(Nfft/2))/Nfft;
        PSD = abs(fft(rx_total,Nfft)/Nfft).^2; PSD = PSD(1:Nfft/2+1);

        plot(ax, f, 10*log10(PSD+eps), 'LineWidth', 1); grid(ax,'on'); hold(ax,'on');
        yl = ylim(ax);

        for k=1:numel(bandsIF)
            xline(ax, bandsIF(k), '-', sprintf('B%d (%.3f GHz)', k, bandsRF(k)/1e9), 'LineWidth',0.6,'Color',[0.4 0.4 0.4]);
        end
        hCF = line(ax, [cf_if cf_if], yl, 'LineStyle','-',  'LineWidth',0.8, 'Color',[0 0.6 0]);
        hJM = line(ax, [jf_if jf_if], yl, 'LineStyle','--', 'LineWidth',0.8, 'Color',[0.85 0 0]);
        try, hCF.Color = [0 0.6 0 0.35]; catch, end
        try, hJM.Color = [0.85 0 0 0.35]; catch, end

        ylim(ax,yl);
        legend(ax,{'PSD','Bands','CF','Jam'},'Location','northeastoutside');
        set([hCF hJM],'HitTest','off','PickableParts','none');
        set(findobj(ax,'Type','ConstantLine'),'HitTest','off','PickableParts','none');
        hold(ax,'off');
        xlabel(ax,'Frekans (IF, Hz)'); ylabel(ax,'PSD'); title(ax,'Spektrum (IF) – RF etiketli bantlar');
    end

    function [llr, rxSym, symbols, modType] = modChanDemod(encBits, EbNo_dB)
        modType = select_modulation(EbNo_dB);
        switch modType
            case 'BPSK'
                symbols = 1 - 2*double(encBits);
                Eb = 1; N0 = Eb / (10^(EbNo_dB/10));
                noise_ch = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
                rxSym = symbols + noise_ch;
                llr = (2/N0)*real(rxSym);
            otherwise
                L = log2(str2double(modType(1:end-3)));
                M = 2^L;
                bits = reshape(encBits, L, []).';
                symbols = qammod(bi2de(bits,'left-msb'), M, 'UnitAveragePower', true);
                Es = mean(abs(symbols).^2); Eb = Es / L;
                N0 = Eb / (10^(EbNo_dB/10));
                noise_ch = sqrt(N0/2)*(randn(size(symbols)) + 1j*randn(size(symbols)));
                rxSym = symbols + noise_ch;
                llr = qamdemod(rxSym, M, 'OutputType','approxllr','UnitAveragePower', true, 'NoiseVariance', N0);
        end
        llr = llr(:);
    end

    function plotConst(ax, tx, rx, modType, EbNo_dB)
        cla(ax); hold(ax,'on');
        n = min(4000, numel(tx));
        plot(ax, real(tx(1:n)), imag(tx(1:n)), '.', 'MarkerSize',6);
        plot(ax, real(rx(1:n)), imag(rx(1:n)), '.', 'MarkerSize',6);
        grid(ax,'on'); xlabel(ax,'I'); ylabel(ax,'Q');
        title(ax, sprintf('Konstelasyon: %s @ %.1f dB',modType,EbNo_dB));
        legend(ax,{'Tx','Rx'},'Location','northeastoutside'); hold(ax,'off');
    end

    function plotLLR(ax, llr)
        cla(ax); histogram(ax, llr, 100); grid(ax,'on');
        xlabel(ax,'LLR'); ylabel(ax,'Adet'); title(ax,'LLR Dağılımı');
    end

    function updateProgress(p), p=max(0,min(1,p)); gauge.Value=p; lblProg.Text=sprintf('%.1f %%',100*p); end
    function logmsg(s), ta.Value=[ta.Value; string(datestr(now,'HH:MM:SS'))+" | "+string(s)]; drawnow limitrate; end

    % helpers (UI metin)
    function txt = fmtBytes(bytes)
        units = {'B','KB','MB','GB','TB'}; v = double(bytes); i=1;
        while v>=1024 && i<numel(units), v=v/1024; i=i+1; end
        txt = sprintf('%.2f %s', v, units{i});
    end
    function txt = fmtTime(sec)
        if ~isfinite(sec) || sec<0, txt='-'; return; end
        sec=round(sec); h=floor(sec/3600); m=floor(mod(sec,3600)/60); s=mod(sec,60);
        if h>0, txt=sprintf('%d saat %d dakika %d saniye',h,m,s);
        elseif m>0, txt=sprintf('%d dakika %d saniye',m,s);
        else, txt=sprintf('%d saniye',s); end
    end
    function y = iff(cond,a,b), if cond, y=a; else, y=b; end, end

    function m = select_modulation(EbNo_dB)
        if EbNo_dB < 4, m='BPSK';
        elseif EbNo_dB < 8, m='4QAM';
        elseif EbNo_dB < 15, m='16QAM';
        else, m='64QAM'; end
    end
end
