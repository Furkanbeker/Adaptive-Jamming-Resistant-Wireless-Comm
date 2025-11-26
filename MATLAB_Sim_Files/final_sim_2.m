function FINAL_SIM()
% MediComm 5-Band Live — Rev8
% İstekler:
% 1) Isı haritası içinde sağ üstte küçük yatay renk çubuğu (Temiz/Jammed)
% 2) Bant tablosu: B1..B5 ve sıkı kolon genişlikleri
% 3) Sağ üst: "İlerleme: %" başlığı + bar + altına Parça/Mod/CF/Jam kutusu
% 4) Spektrumda CF/Jam çizgileri kalın

%% --------- ŞARTNAME SABİTLERİ ---------
BANDS_RF_HZ = [2.402e9, 2.406e9, 2.410e9, 2.414e9, 2.418e9];   % RF merkezler
BANDS_IF_HZ = [100e3,   200e3,   300e3,   400e3,   500e3];     % IF merkezler

%% --------- STATE ----------
state.running       = false;
state.paused        = false;
state.stopRequested = false;
state.simInit       = false;

state.videoPath     = "";
state.outputFolder  = "";

state.fs            = 1.2e6;       % IF örnekleme
state.T             = 1;           % görsel kanal penceresi
state.bandsRF       = BANDS_RF_HZ(:).';
state.bandsIF       = BANDS_IF_HZ(:).';

state.jammer_amp    = 2.0;
state.EbNo_dB       = 6;
state.snrAuto       = false;

state.chunk_size    = 100000;      % B
state.maxItersLDPC  = 50;

state.jamAutoWander = false;
state.jamIdx        = 3;           % jammer bandı
state.cfIdx         = 3;           % çalışma bandı

state.darkTheme     = false;

state.all_rx_bytes  = uint8([]);
state.csvEnable     = true;
state.csvRows       = {};

% Çalışma oturumu (resume)
state.sim = struct('video_bytes',[], 'total_bytes',0, 'num_chunks',0, ...
                   'outFile',"", 'part_idx',0, 'elapsedActive',0, ...
                   'berPreSeries',[], 'berPostSeries',[], 'heat',[]);

% LDPC/CRC
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
state.K = encCfg.NumInformationBits;
state.N = encCfg.BlockLength;
state.crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
state.crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');

%% --------- UI ----------
fig = uifigure('Name','MediComm 5-Band Live','Position',[80 60 1600 900]);
fig.WindowKeyPressFcn = @(~,ev) onKey(ev);

gl = uigridlayout(fig,[1,3], 'ColumnWidth',{360,'1x',420}, 'RowHeight',{'1x'}, ...
    'Padding',8,'ColumnSpacing',10,'RowSpacing',8);

%% --- SOL: Kontroller ---
pLeft = uipanel(gl,'Title','Kontroller','FontWeight','bold'); pLeft.Layout.Column=1;
gLeft = uigridlayout(pLeft,[3,1],'RowHeight',{'fit','fit','fit'},'Padding',6,'RowSpacing',8);

% 1) Dosya & Çalıştır
pIO = uipanel(gLeft,'Title','Dosya & Çalıştır'); pIO.Layout.Row=1;
gIO = uigridlayout(pIO,[3,3], 'ColumnWidth',{110, 110, 110}, 'RowHeight',{'fit','fit','fit'}, 'Padding',6, 'ColumnSpacing',8,'RowSpacing',6);
uibutton(gIO,'Text','Video Seç','ButtonPushedFcn',@(b,~) onChooseVideo());
lblVideo = uilabel(gIO,'Text','yok','WordWrap','on'); lblVideo.Layout.Column=[2 3];
uibutton(gIO,'Text','Çıkış klasörü','ButtonPushedFcn',@(b,~) onChooseOutput());
lblOut = uilabel(gIO,'Text','(varsayılan: girişle aynı)','WordWrap','on'); lblOut.Layout.Column=[2 3];
btnStart = uibutton(gIO,'Text','BAŞLAT','BackgroundColor',[0.21 0.6 0.24], ...
    'FontWeight','bold','ButtonPushedFcn',@(b,~) onStartNew());
btnStop  = uibutton(gIO,'Text','DURDUR','BackgroundColor',[0.75 0.2 0.2], ...
    'Enable','off','ButtonPushedFcn',@(b,~) onStopOrResume());
cbTheme = uicheckbox(gIO,'Text','Koyu Tema','Value',state.darkTheme,'ValueChangedFcn',@(c,~) toggleTheme(c.Value));

% 2) Kanal & SNR
pSNR = uipanel(gLeft,'Title','Kanal & SNR'); pSNR.Layout.Row=2;
gS = uigridlayout(pSNR,[3,2],'ColumnWidth',{90,'1x'},'RowHeight',{'fit','fit','fit'},'Padding',6,'RowSpacing',6,'ColumnSpacing',10);
uilabel(gS,'Text','SNR (dB):','HorizontalAlignment','right');
sSNR = uislider(gS,'Limits',[-2 20],'Value',state.EbNo_dB,'MajorTicks',-2:2:20, ...
    'MinorTicks',[],'ValueChangedFcn',@(s,~) setSNR(s.Value));
uilabel(gS,'Text','SNR Oto:','HorizontalAlignment','right');
cbSNRAuto = uicheckbox(gS,'Text','Rastgele (parça başı)','Value',state.snrAuto,'ValueChangedFcn',@(c,~) setSNRAuto(c.Value));
uibutton(gS,'Text','SNR rastgele','ButtonPushedFcn',@(b,~) randomizeSNRonce()); 

% 3) Jammer
pJam = uipanel(gLeft,'Title','Jammer'); pJam.Layout.Row=3;
gJ = uigridlayout(pJam,[3,2],'ColumnWidth',{90,'1x'},'RowHeight',{'fit','fit','fit'},'Padding',6,'RowSpacing',6,'ColumnSpacing',10);
uilabel(gJ,'Text','Güç (amp):','HorizontalAlignment','right');
sJamA = uislider(gJ,'Limits',[0 5],'Value',state.jammer_amp,'MajorTicks',0:1:5,'MinorTicks',[], ...
    'ValueChangedFcn',@(s,~) setJamAmp(s.Value));

uilabel(gJ,'Text','Band:','HorizontalAlignment','right');
jamBtns = uigridlayout(gJ,[5,1], 'RowHeight',{36,36,36,36,36}, 'ColumnWidth',{'1x'}, 'RowSpacing',6, 'Padding',4);
btnJ = gobjects(1,5);
for k = 1:5
    txt = sprintf('B%d (%.3f)', k, state.bandsRF(k)/1e9);
    btnJ(k) = uibutton(jamBtns,'Text',txt,'ButtonPushedFcn', @(~,~) selectJammer(k));
    btnJ(k).Layout.Row = k; btnJ(k).Layout.Column = 1;
end
uilabel(gJ,'Text','Jammer Oto:','HorizontalAlignment','right');
cbWander = uicheckbox(gJ,'Text','Band gezer','Value',state.jamAutoWander, 'ValueChangedFcn',@(c,~) setWander(c.Value));

%% --- ORTA: GRAFİKLER ---
pMid = uipanel(gl,'Title','Grafikler','FontWeight','bold'); pMid.Layout.Column=2;
gMid = uigridlayout(pMid,[2,2],'RowHeight',{'1x','1x'},'ColumnWidth',{'1x','1x'},'Padding',6,'RowSpacing',8,'ColumnSpacing',10);

axSpec = uiaxes(gMid); axSpec.Layout.Row=1; axSpec.Layout.Column=1;
title(axSpec,'Spektrum (IF: tıkla snap, çift tık büyüt)'); xlabel(axSpec,'Frekans (IF, Hz)'); ylabel(axSpec,'PSD');
axSpec.ButtonDownFcn = @(~,~) onSpectrumClick();

axBER = uiaxes(gMid); axBER.Layout.Row=1; axBER.Layout.Column=2; hold(axBER,'on');
title(axBER,'BER Zaman Serisi'); xlabel(axBER,'Parça #'); ylabel(axBER,'BER');
pltPre  = plot(axBER,nan,nan,'-'); 
pltPost = plot(axBER,nan,nan,'-'); 
legend(axBER,{'Pre-LDPC','Post-LDPC'},'Location','northeast');

axConst = uiaxes(gMid); axConst.Layout.Row=2; axConst.Layout.Column=1;
title(axConst,'Konstelasyon'); xlabel(axConst,'I'); ylabel(axConst,'Q');

axHeat = uiaxes(gMid); axHeat.Layout.Row=2; axHeat.Layout.Column=2;
title(axHeat,'Bant Isı Haritası'); xlabel(axHeat,'Parça #'); ylabel(axHeat,'Band');
colormap(axHeat, [0 0.6 0; 0.85 0 0]); % 0=Temiz, 1=Jammed
% kendi küçük renk çubuğunu çizmek için ilk kurulum
drawHeatLegend(axHeat,true);  % boş olsa da yer tutucu

%% --- SAĞ: METRİKLER & KAYIT ---
pRight = uipanel(gl,'Title','Metrikler & Kayıt','FontWeight','bold'); pRight.Layout.Column=3;
% 8 satır: Durum | "İlerleme: %" | Bar | Parça/Mod/CF/Jam kutusu | Xfer | Tablo | Log | Kısayol
gRight = uigridlayout(pRight,[8,1], 'RowHeight',{'fit','fit','fit','fit','1x','1x','1x','fit'}, 'Padding',6,'RowSpacing',8);

% 0) durum
lblStats = uilabel(gRight,'Text','hazır','WordWrap','on');

% 1) başlık
lblProgTitle = uilabel(gRight,'Text','İlerleme: 0.0 %','FontWeight','bold');

% 2) progress bar
pbOuter = uipanel(gRight,'BackgroundColor',[0.95 0.95 0.95],'BorderType','line');
pbInner = uipanel(pbOuter,'BackgroundColor',[0.2 0.6 0.9],'Position',[3 3 0 16]);

% 3) anlık durum kutusu
pNow = uipanel(gRight,'Title','Anlık Durum'); 
gNow = uigridlayout(pNow,[4,2],'ColumnWidth',{70,'1x'},'RowHeight',{'fit','fit','fit','fit'},'Padding',6,'RowSpacing',4,'ColumnSpacing',8);
uilabel(gNow,'Text','Parça:','HorizontalAlignment','right');  lblNowPart = uilabel(gNow,'Text','-');
uilabel(gNow,'Text','Mod:','HorizontalAlignment','right');    lblNowMod  = uilabel(gNow,'Text','-');
uilabel(gNow,'Text','CF:','HorizontalAlignment','right');     lblNowCF   = uilabel(gNow,'Text','-');
uilabel(gNow,'Text','Jam:','HorizontalAlignment','right');    lblNowJam  = uilabel(gNow,'Text','-');

% 4) aktarım durumu
pXfer = uipanel(gRight,'Title','Aktarım Durumu'); 
gx = uigridlayout(pXfer,[5,2],'ColumnWidth',{120,'1x'},'RowHeight',repmat({'fit'},1,5), ...
    'ColumnSpacing',8,'RowSpacing',4,'Padding',6);
uilabel(gx,'Text','Dosya Boyutu:','HorizontalAlignment','right');  lblSize   = uilabel(gx,'Text','-');
uilabel(gx,'Text','İletilen:','HorizontalAlignment','right');      lblSent   = uilabel(gx,'Text','-');
uilabel(gx,'Text','Kalan:','HorizontalAlignment','right');         lblRemain = uilabel(gx,'Text','-');
uilabel(gx,'Text','Hız:','HorizontalAlignment','right');           lblRate   = uilabel(gx,'Text','-');
uilabel(gx,'Text','Kalan Süre:','HorizontalAlignment','right');    lblETA    = uilabel(gx,'Text','-');

% 5) 5 Bant Taraması Tablosu — dar kolonlar ve B1..B5
tbl = uitable(gRight,'ColumnName',{'Band','RF_GHz','Power_dB','BG_dB','Thresh','Durum'},'Data',cell(0,6));
tbl.ColumnWidth = {45,60,70,68,64,60};  % sıkılaştır
tbl.RowStriping = 'on';

% 6) Log & CSV & Mod
pLog = uipanel(gRight,'Title','Log & Kaydet');
gLog = uigridlayout(pLog,[1,2],'ColumnWidth',{'1x',120},'Padding',6,'ColumnSpacing',8);
ta = uitextarea(gLog,'Editable','off','FontName','Consolas');
pBadge = uipanel(gLog,'BorderType','line','Title','Mod'); 
gBadge = uigridlayout(pBadge,[2,1],'RowHeight',{'fit','fit'},'Padding',6,'RowSpacing',6);
lblMod = uilabel(gBadge,'Text','  -  ','FontWeight','bold','HorizontalAlignment','center','BackgroundColor',[0.92 0.92 0.92]);
cbCSV = uicheckbox(gBadge,'Text','CSV log kaydet','Value',state.csvEnable,'ValueChangedFcn',@(c,~) setCSV(c.Value));

% 7) alt yardım
lblHelp = uilabel(gRight,'Text','Kısayollar: Space=Duraklat/Devam | J=Jammer sırayla | S=SNR rastgele (tek) | D=Koyu Tema','HorizontalAlignment','center');

applyTheme();

%% --------- CALLBACKS (UI) ----------
    function onKey(ev)
        switch lower(ev.Key)
            case 'space', onStopOrResume();
            case 'j', selectJammer(mod(state.jamIdx,5)+1);
            case 's', randomizeSNRonce();
            case 'd', toggleTheme(~state.darkTheme); cbTheme.Value = state.darkTheme;
        end
    end

    function onSpectrumClick()
        if strcmp(fig.SelectionType,'open') % çift tık: büyük pencerede
            popSpectrum(); return;
        end
        cp = axSpec.CurrentPoint; x = cp(1,1);
        [~,idx] = min(abs(state.bandsIF - x));
        state.cfIdx = idx; logmsg(sprintf("CF grafikten: B%d (%.3f GHz)",idx,state.bandsRF(idx)/1e9));
    end

    function popSpectrum()
        f2 = figure('Name','Spektrum (Büyük)','Color','w');
        fs=state.fs; t = 0:1/fs:state.T-1/fs;
        cf=state.bandsIF(state.cfIdx); jf=state.bandsIF(state.jamIdx);
        rx = sin(2*pi*cf*t) + state.jammer_amp*sin(2*pi*jf*t) + randn(size(t))*0.1;
        Nfft = 2^14; f = fs*(0:(Nfft/2))/Nfft; P = abs(fft(rx,Nfft)/Nfft).^2; P=P(1:Nfft/2+1);
        plot(f,10*log10(P+eps),'LineWidth',1); grid on; hold on; yl=ylim;
        plot([cf cf],yl,'-','LineWidth',2.5,'Color',[0 0.6 0]);
        plot([jf jf],yl,'--','LineWidth',2.5,'Color',[0.85 0 0]); hold off;
        xlabel('Frekans (IF, Hz)'); ylabel('PSD'); title('Spektrum (Büyük)');
    end

    function onChooseVideo()
        [f,p] = uigetfile({'*.mp4;*.bin;*.dat','Video/Raw (*.mp4,*.bin,*.dat)';'*.*','Hepsi'});
        if isequal(f,0), return; end
        state.videoPath = fullfile(p,f); lblVideo.Text = f;
        logmsg("Video seçildi: "+state.videoPath);
    end

    function onChooseOutput()
        p = uigetdir();
        if isequal(p,0), return; end
        state.outputFolder = p; lblOut.Text = p; logmsg("Çıkış klasörü: "+p);
    end

    function setSNR(v), state.EbNo_dB=v; end
    function setSNRAuto(tf), state.snrAuto=tf; sSNR.Enable=tern(tf,'off','on'); end
    function setJamAmp(v), state.jammer_amp=v; end
    function setWander(tf), state.jamAutoWander=tf; for k=1:5, btnJ(k).Enable=tern(tf,'off','on'); end, end
    function setCSV(tf), state.csvEnable=tf; end

    function selectJammer(k)
        state.jamIdx=k; logmsg(sprintf("Jammer band: B%d (%.3f GHz)",k,state.bandsRF(k)/1e9));
    end
    function randomizeSNRonce()
        v = round(sSNR.Limits(1) + (sSNR.Limits(2)-sSNR.Limits(1))*rand(),2);
        state.EbNo_dB=v; sSNR.Value=v; logmsg(sprintf("SNR tek sefer: %.2f dB",v));
    end
    function toggleTheme(tf), state.darkTheme=tf; applyTheme(); end

    function onStartNew()
        if state.videoPath=="", uialert(fig,'Lütfen video dosyası seçin.','Uyarı'); return; end
        if state.running, return; end
        initSimulation();
        state.paused=false; state.stopRequested=false; state.running=true;
        btnStart.Enable='off'; btnStop.Enable='on'; btnStop.Text='DURDUR'; btnStop.BackgroundColor=[0.75 0.2 0.2];
        lblStats.Text='çalışıyor...';
        loopSimulation(state.sim.part_idx+1);
        state.running=false;
        if state.paused, lblStats.Text='durakladı';
        else, lblStats.Text='bitti'; btnStart.Enable='on'; btnStop.Enable='off'; btnStop.Text='DURDUR'; end
    end

    function onStopOrResume()
        if state.running % DURDUR
            state.stopRequested=true; state.paused=true;
            btnStop.Text='DEVAM ET'; btnStop.BackgroundColor=[0.6 0.6 0.6]; lblStats.Text='duraklatılıyor...';
        else % DEVAM
            if ~state.simInit || ~state.paused, return; end
            state.stopRequested=false; state.paused=false; state.running=true;
            btnStop.Text='DURDUR'; btnStop.BackgroundColor=[0.75 0.2 0.2]; btnStart.Enable='off'; lblStats.Text='çalışıyor...';
            loopSimulation(state.sim.part_idx+1);
            state.running=false;
            if state.paused, lblStats.Text='durakladı';
            else, lblStats.Text='bitti'; btnStart.Enable='on'; btnStop.Enable='off'; btnStop.Text='DURDUR'; end
        end
    end

%% --------- ÇEKİRDEK ----------
    function initSimulation()
        [video_bytes,total_bytes] = readBytes(state.videoPath);
        state.all_rx_bytes = uint8([]);
        state.sim.video_bytes = video_bytes;
        state.sim.total_bytes = total_bytes;
        state.sim.num_chunks  = ceil(double(total_bytes)/double(state.chunk_size));
        state.sim.outFile     = resolveOutFile();
        state.sim.part_idx    = 0; state.sim.elapsedActive=0;
        state.sim.berPreSeries= nan(1,state.sim.num_chunks);
        state.sim.berPostSeries=nan(1,state.sim.num_chunks);
        state.sim.heat        = nan(5,state.sim.num_chunks);
        state.simInit = true;

        % reset metrikler
        lblSize.Text=fmtBytes(total_bytes); lblSent.Text=fmtBytes(0); lblRemain.Text=fmtBytes(total_bytes);
        lblRate.Text='-'; lblETA.Text='-'; updateProgress(0);
        set(pltPre,'XData',nan,'YData',nan); set(pltPost,'XData',nan,'YData',nan);
        cla(axHeat); title(axHeat,'Bant Isı Haritası'); ylabel(axHeat,'Band'); xlabel(axHeat,'Parça #');
        drawHeatLegend(axHeat,true);

        % CSV header
        state.csvRows = [{'timestamp','part','snr_db','jam_rf_GHz','cf_rf_GHz','ber_pre','ber_post','rate_Mbps','eta_s'}];
    end

    function loopSimulation(start_part)
        for part = start_part : state.sim.num_chunks
            if state.stopRequested, break; end
            ticStep = tic;

            % SNR otomatik (parça başı)
            if state.snrAuto
                state.EbNo_dB = round(sSNR.Limits(1)+(sSNR.Limits(2)-sSNR.Limits(1))*rand(),2);
                sSNR.Value=state.EbNo_dB;
            end
            EbNo_dB = state.EbNo_dB; jamAmp = state.jammer_amp;

            % Jammer otomatik gezer
            if state.jamAutoWander, state.jamIdx = randi([1 5]); end
            jf_if = state.bandsIF(state.jamIdx);
            cf_if = state.bandsIF(state.cfIdx);

            % Parça sınırları
            idx_start = (part-1)*state.chunk_size + 1;
            idx_end   = min(part*state.chunk_size, state.sim.total_bytes);
            part_bytes= state.sim.video_bytes(idx_start:idx_end);

            % Encode + kanal + decode (özet)
            video_bits = reshape(de2bi(part_bytes,8,'left-msb')',[],1);
            txCRC = state.crcGen(logical(video_bits));
            pad = mod(state.K - mod(length(txCRC),state.K), state.K);
            if pad>0, txCRC = [txCRC; false(pad,1)]; end
            numBl = length(txCRC)/state.K;
            encBits = false(numBl*state.N,1);
            for i=1:numBl, encBits((i-1)*state.N+1:i*state.N)=ldpcEncode(txCRC((i-1)*state.K+1:i*state.K),encCfg); end

            % 5 bant taraması + kaçış
            [clean_if, tableRows, jamVec] = scanJammer5(jf_if, jamAmp, state.bandsIF, state.bandsRF, state.fs, state.T, EbNo_dB);
            tbl.Data = tableRows;                       % tabloyu göster
            state.sim.heat(:,part) = jamVec(:);         % ısı haritası
            if ~isempty(clean_if)
                [~,idx] = min(abs(state.bandsIF-clean_if(1))); state.cfIdx = idx; cf_if = state.bandsIF(idx);
            end

            % Görsel kanal -> Spektrum
            t = 0:1/state.fs:state.T-1/state.fs;
            rx_total = sin(2*pi*cf_if*t) + jamAmp*sin(2*pi*jf_if*t) + sqrt(1/(10^(EbNo_dB/10))/2)*randn(size(t));
            plotPSD(axSpec, rx_total, state.fs, state.bandsIF, state.bandsRF, cf_if, jf_if);

            % Mod/Demod + BER
            [llr, rxSym, symbols, modType] = modChanDemod(encBits, EbNo_dB);
            plotConst(axConst, symbols, rxSym, modType, EbNo_dB);

            rxHard = llr(:) < 0; ber_pre = mean(rxHard ~= encBits);
            decBits = false(numBl*state.K,1);
            for i=1:numBl, decBits((i-1)*state.K+1:i*state.K) = ldpcDecode(llr((i-1)*state.N+1:i*state.N), decCfg, state.maxItersLDPC); end
            ber_post = mean(decBits(1:length(txCRC)) ~= txCRC);

            % çıktı biriktir (sim: orijinal parça)
            state.all_rx_bytes = [state.all_rx_bytes; uint8(part_bytes)]; %#ok<AGROW>

            % aktarım metrikleri
            state.sim.part_idx = part;
            state.sim.elapsedActive = state.sim.elapsedActive + toc(ticStep);
            sent_bytes   = double(idx_end);
            remain_bytes = double(state.sim.total_bytes - idx_end);
            rate_Bps = sent_bytes / max(state.sim.elapsedActive, eps);
            rate_Mbps = (rate_Bps*8)/1e6;
            eta_s  = remain_bytes / max(rate_Bps, eps);

            % sağ panel metinleri
            lblSize.Text=fmtBytes(state.sim.total_bytes);
            lblSent.Text=fmtBytes(sent_bytes);
            lblRemain.Text=fmtBytes(remain_bytes);
            lblRate.Text=sprintf('%.2f Mbit/s (%.2f MB/s)', rate_Mbps, rate_Bps/1e6);
            lblETA.Text = fmtTime(eta_s);
            lblNowPart.Text = sprintf('%d / %d', part, state.sim.num_chunks);
            lblNowMod.Text  = modType;
            lblNowCF.Text   = sprintf('%.3f GHz', state.bandsRF(state.cfIdx)/1e9);
            lblNowJam.Text  = sprintf('%.3f GHz', state.bandsRF(state.jamIdx)/1e9);

            % BER serileri & ilerleme
            state.sim.berPreSeries(part)=ber_pre; state.sim.berPostSeries(part)=ber_post;
            set(pltPre,'XData',1:part,'YData',state.sim.berPreSeries(1:part));
            set(pltPost,'XData',1:part,'YData',state.sim.berPostSeries(1:part));
            updateProgress(double(idx_end)/double(state.sim.total_bytes));

            % ısı haritasını güncelle
            imagesc(axHeat, 1:part, 1:5, state.sim.heat(:,1:part)); axis(axHeat,'xy');
            yticks(axHeat,1:5); yticklabels(axHeat,{'B1','B2','B3','B4','B5'});
            drawHeatLegend(axHeat,false);  % ufak yatay çubuk

            % CSV
            if state.csvEnable
                state.csvRows(end+1,:) = {datestr(now,'HH:MM:SS.FFF'), part, EbNo_dB, ...
                    state.bandsRF(state.jamIdx)/1e9, state.bandsRF(state.cfIdx)/1e9, ...
                    ber_pre, ber_post, rate_Mbps, eta_s}; %#ok<AGROW>
            end

            lblStats.Text = sprintf('Parça %d/%d | Mod: %s | CF=%.3f GHz | Jam=%.3f GHz', ...
                part, state.sim.num_chunks, modType, state.bandsRF(state.cfIdx)/1e9, state.bandsRF(state.jamIdx)/1e9);
            drawnow limitrate;

            if state.stopRequested, break; end
        end

        % tamamlandıysa dosya + csv yaz
        if state.sim.part_idx >= state.sim.num_chunks
            writeOutputAndCSV();
        end
    end

%% --------- HELPERLAR ----------
    function writeOutputAndCSV()
        if ~isempty(state.all_rx_bytes)
            try
                fid=fopen(state.sim.outFile,'wb'); fwrite(fid,state.all_rx_bytes,'uint8'); fclose(fid);
                logmsg("Çıktı yazıldı: "+state.sim.outFile);
            catch ME, logmsg("Çıkış yazım hatası: "+ME.message);
            end
        end
        if state.csvEnable && numel(state.csvRows)>1
            [outDir,~,~]=fileparts(state.sim.outFile);
            csvPath=fullfile(outDir, ['session_', datestr(now,'yyyymmdd_HHMMSS'), '.csv']);
            try
                fh=fopen(csvPath,'w'); fprintf(fh,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n', state.csvRows{1,:});
                for r=2:size(state.csvRows,1)
                    row=state.csvRows(r,:); fprintf(fh,'%s,%d,%.2f,%.3f,%.3f,%.6f,%.6f,%.2f,%.0f\n', row{:});
                end
                fclose(fh); logmsg("CSV yazıldı: "+csvPath);
            catch ME, logmsg("CSV yazım hatası: "+ME.message);
            end
        end
    end

    function outFile = resolveOutFile()
        if state.outputFolder=="", [p,~,~]=fileparts(state.videoPath); outFolder=p; else, outFolder=state.outputFolder; end
        outFile = fullfile(outFolder, ['received_video_', datestr(now,'yyyymmdd_HHMMSS'), '.mp4']);
    end

    function [video_bytes,total_bytes] = readBytes(path)
        try, fid=fopen(path,'rb'); if fid==-1, error('Dosya açılamadı'); end
            video_bytes=fread(fid,'uint8'); fclose(fid);
            total_bytes=length(video_bytes); logmsg("Video boyutu: "+fmtBytes(total_bytes));
        catch
            logmsg("Dosya okuma hatası. Sentetik veri (1 MB) oluşturuluyor...");
            video_bytes=uint8(randi([0 255],1e6,1)); total_bytes=numel(video_bytes);
        end
    end

    % ----- 5 bant taraması; tablo "B1..B5" ve dar kolonlar -----
    function [clean_if, tableRows, jamVec] = scanJammer5(jf_if, jam_amp, bandsIF, bandsRF, fs, T, EbNo_dB)
        t=0:1/fs:T-1/fs; N0=1/(10^(EbNo_dB/10));
        noise=sqrt(N0/2)*randn(size(t));
        tableRows=cell(0,6); clean_if=[]; jamVec=zeros(5,1);
        rxJam = jam_amp*sin(2*pi*jf_if*t)+noise;
        Nfft=2^nextpow2(length(rxJam)); f=fs*(0:(Nfft/2))/Nfft;
        PSD=abs(fft(rxJam,Nfft)/Nfft).^2; PSD=PSD(1:Nfft/2+1);
        bg_band=PSD(1:round(numel(PSD)*0.25)); bg_dB=10*log10(mean(bg_band)+eps);
        for k=1:5
            cf=bandsIF(k);
            idx=abs(f-cf)<1000; % 1 kHz pencere
            if ~any(idx), pwr_dB=-Inf; else
                seg=PSD(idx); seg=sort(seg,'descend'); pwr_dB=10*log10(mean(seg(1:min(5,numel(seg))))+eps);
            end
            thr=bg_dB+15;
            if pwr_dB>thr, stat='JAMMED'; jamVec(k)=1; else, stat='TEMIZ'; clean_if(end+1)=cf; end %#ok<AGROW>
            % B1..B5 ve kısa RF (GHz)
            tableRows(end+1,:)={sprintf('B%d',k), sprintf('%.3f',bandsRF(k)/1e9), pwr_dB, bg_dB, thr, stat}; %#ok<AGROW>
        end
    end

    % ----- Spektrum (CF/Jam kalın) -----
    function plotPSD(ax, rx, fs, bandsIF, bandsRF, cf, jf)
        Nfft=2^nextpow2(length(rx)); f=fs*(0:(Nfft/2))/Nfft;
        P=abs(fft(rx,Nfft)/Nfft).^2; P=P(1:Nfft/2+1);
        plot(ax, f, 10*log10(P+eps), 'LineWidth',1.0); grid(ax,'on'); hold(ax,'on');
        yl=ylim(ax);
        for k=1:numel(bandsIF)
            xline(ax,bandsIF(k),'-',sprintf('B%d (%.3f GHz)',k,bandsRF(k)/1e9),'LineWidth',0.7);
        end
        plot(ax,[cf cf],yl,'-','LineWidth',2.5,'Color',[0 0.6 0]);    % daha kalın CF
        plot(ax,[jf jf],yl,'--','LineWidth',2.5,'Color',[0.85 0 0]);  % daha kalın Jam
        ylim(ax,yl); hold(ax,'off'); xlabel(ax,'Frekans (IF, Hz)'); ylabel(ax,'PSD'); title(ax,'Spektrum (IF)');
    end

    % ----- Mod/Demod -----
    function [llr, rxSym, symbols, modType] = modChanDemod(encBits, EbNo_dB)
        modType = select_modulation(EbNo_dB);
        switch modType
            case 'BPSK'
                symbols = 1 - 2*double(encBits);
                Eb = 1; N0 = Eb / (10^(EbNo_dB/10));
                noise_ch = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
                rxSym = symbols + noise_ch; llr = (2/N0)*real(rxSym);
            otherwise
                L = log2(str2double(modType(1:end-3))); M=2^L;
                bits = reshape(encBits,L,[]).';
                symbols = qammod(bi2de(bits,'left-msb'),M,'UnitAveragePower',true);
                Es = mean(abs(symbols).^2); Eb=Es/L; N0=Eb/(10^(EbNo_dB/10));
                noise_ch = sqrt(N0/2)*(randn(size(symbols))+1j*randn(size(symbols)));
                rxSym = symbols + noise_ch;
                llr = qamdemod(rxSym,M,'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',N0);
        end
        llr=llr(:);
    end

    function plotConst(ax, tx, rx, modType, EbNo)
        cla(ax); hold(ax,'on');
        n=min(4000,numel(tx));
        plot(ax,real(tx(1:n)),imag(tx(1:n)),'.','MarkerSize',6);
        plot(ax,real(rx(1:n)),imag(rx(1:n)),'.','MarkerSize',6);
        grid(ax,'on'); xlabel(ax,'I'); ylabel(ax,'Q');
        title(ax,sprintf('Konstelasyon: %s @ %.1f dB',modType,EbNo));
        legend(ax,{'Tx','Rx'},'Location','northeast'); 
        setBadge(modType);
    end

    function setBadge(modType)
        lblMod.Text = ['  ' modType '  '];
        switch modType
            case 'BPSK',  lblMod.BackgroundColor=[0.60 0.86 0.60];
            case '4QAM',  lblMod.BackgroundColor=[0.60 0.76 0.95];
            case '16QAM', lblMod.BackgroundColor=[0.98 0.84 0.60];
            case '64QAM', lblMod.BackgroundColor=[0.95 0.60 0.60];
        end
        lblNowMod.Text = modType;
    end

    function updateProgress(p)
        p=max(0,min(1,p));
        lblProgTitle.Text = sprintf('İlerleme: %.1f %%',100*p);
        outerW = pbOuter.Position(3); pbInner.Position=[3 3 max(0,(outerW-6)*p) 16];
    end

    function applyTheme()
        if state.darkTheme
            fig.Color=[0.10 0.11 0.12];
            for c=[pLeft pMid pRight], c.BackgroundColor=[0.13 0.14 0.15]; end
            for ax=[axSpec axBER axConst axHeat], ax.Color=[0.17 0.18 0.2]; ax.GridColor=[0.6 0.6 0.6]; ax.XColor=[0.95 0.95 0.95]; ax.YColor=[0.95 0.95 0.95]; end
        else
            fig.Color='w';
            for c=[pLeft pMid pRight], c.BackgroundColor='w'; end
            for ax=[axSpec axBER axConst axHeat], ax.Color='w'; ax.GridColor=[0.35 0.35 0.35]; ax.XColor='k'; ax.YColor='k'; end
        end
    end

    function logmsg(s)
        ta.Value = [ta.Value; string(datestr(now,'HH:MM:SS')) + " | " + string(s)];
        drawnow limitrate;
    end

    % helpers (format)
    function txt = fmtBytes(bytes)
        units={'B','KB','MB','GB','TB'}; v=double(bytes); i=1; while v>=1024 && i<numel(units), v=v/1024; i=i+1; end
        txt=sprintf('%.2f %s',v,units{i});
    end
    function txt = fmtTime(sec)
        if ~isfinite(sec) || sec<0, txt='-'; return; end
        sec=round(sec); h=floor(sec/3600); m=floor(mod(sec,3600)/60); s=mod(sec,60);
        if h>0, txt=sprintf('%d saat %d dakika %d saniye',h,m,s);
        elseif m>0, txt=sprintf('%d dakika %d saniye',m,s);
        else, txt=sprintf('%d saniye',s); end
    end
    function y=tern(c,a,b), if c, y=a; else, y=b; end, end
    function m=select_modulation(EbNo)
        if EbNo<4, m='BPSK'; elseif EbNo<8, m='4QAM'; elseif EbNo<15, m='16QAM'; else, m='64QAM'; end
    end

    % ----- Isı haritası içinde sağ üstte ufak yatay renk çubuğu -----
    function drawHeatLegend(ax,first)
        % önceki legend öğelerini sil
        delete(findobj(ax,'Tag','HeatLegend'));
        % eksen limitleri
        xl = xlim(ax); yl = ylim(ax);
        xs = xl(2)-xl(1); ys = yl(2)-yl(1);
        % küçük kutu boyutu ve konumu (sağ üst)
        w = 0.28*xs; h = 0.16*ys;
        x0 = xl(2) - w - 0.03*xs;
        y0 = yl(2) - h - 0.03*ys;
        % arka plan
        rectangle(ax,'Position',[x0 y0 w h],'FaceColor',[1 1 1 0.75],'EdgeColor',[0.3 0.3 0.3],'Tag','HeatLegend');
        % sol (Temiz) ve sağ (Jammed) yarılar
        rectangle(ax,'Position',[x0 y0 w/2 h*0.55],'FaceColor',[0 0.6 0],'EdgeColor','none','Tag','HeatLegend');
        rectangle(ax,'Position',[x0+w/2 y0 w/2 h*0.55],'FaceColor',[0.85 0 0],'EdgeColor','none','Tag','HeatLegend');
        % metinler
        text(ax,x0+w*0.25,y0+h*0.75,'Temiz','HorizontalAlignment','center','FontSize',10,'Tag','HeatLegend');
        text(ax,x0+w*0.75,y0+h*0.75,'Jammed','HorizontalAlignment','center','FontSize',10,'Tag','HeatLegend');
        if first, drawnow; end
    end
end
