function ldpc_jammer_gui_allinone()
% MediComm 5-Band Live (RF=2.4 GHz, IF simülasyon)
% - Çalışma bandı slider YOK
% - Jammer 5 buton (B1..B5) ile seçiliyor
% - RF bandları: 2402/2406/2410/2414/2418 MHz
% - İç simülasyon IF'te (100/200/300/400/500 kHz) yürür; UI RF gösterir

%% -------------------- SABİTLER (ŞARTNAME) --------------------
% RF (gerçek) merkezler [Hz]
BANDS_RF_HZ = [2.402e9, 2.406e9, 2.410e9, 2.414e9, 2.418e9];
% IF (simülasyon) merkezleri [Hz] (görsel/algılama burada)
BANDS_IF_HZ = [100e3, 200e3, 300e3, 400e3, 500e3];

%% -------------------- STATE --------------------
state.running       = false;
state.stopRequested = false;
state.videoPath     = "";
state.outputFolder  = "";
state.fs            = 1e6;       % IF örnekleme
state.T             = 1;         % 1 sn pencere
state.bandsRF       = BANDS_RF_HZ(:).';
state.bandsIF       = BANDS_IF_HZ(:).';
state.jammer_amp    = 2.0;
state.EbNo_dB       = 6;
state.chunk_size    = 100000;
state.maxItersLDPC  = 50;
state.progress      = 0;
state.all_rx_bytes  = uint8([]);
state.jamAutoWander = false;     % buton istendiği için default kapalı
state.jamIdx        = 3;         % jammer başlangıç bandı
state.cfIdx         = 3;         % çalışma bandı (kaçışla değişir)

% LDPC/CRC
H = dvbs2ldpc(1/2);
encCfg = ldpcEncoderConfig(H);
decCfg = ldpcDecoderConfig(H);
state.K = encCfg.NumInformationBits;
state.N = encCfg.BlockLength;
state.crcGen = comm.CRCGenerator('Polynomial','0x04C11DB7');
state.crcDet = comm.CRCDetector('Polynomial','0x04C11DB7');

%% -------------------- UI --------------------
fig = uifigure('Name','MediComm 5-Band Live','Position',[60 60 1400 850]);
gl  = uigridlayout(fig,[3,3], 'RowHeight',{'fit','1x','fit'}, ...
    'ColumnWidth',{360,'1x',500}, 'Padding',8,'RowSpacing',8,'ColumnSpacing',10);

% Left controls
pLeft = uipanel(gl,'Title','Kontroller','FontWeight','bold');
pLeft.Layout.Row = [1 3]; pLeft.Layout.Column = 1;
gL = uigridlayout(pLeft,[16,2], 'RowHeight',repmat({'fit'},1,16), ...
    'ColumnWidth',{'1x','1x'}, 'RowSpacing',8,'ColumnSpacing',8,'Padding',8);

uibutton(gL,'Text','Video Seç (mp4/bin)','ButtonPushedFcn',@(btn,~) onChooseVideo());
uilabel(gL,'Text',' ');

uibutton(gL,'Text','Çıkış klasörü','ButtonPushedFcn',@(btn,~) onChooseOutput());
lblOut = uilabel(gL,'Text','Çıkış: (varsayılan: girişle aynı klasör)','WordWrap','on');
lblOut.Layout.Column = [1 2];

uilabel(gL,'Text','SNR (dB)');
sSNR = uislider(gL,'Limits',[-2 20],'Value',state.EbNo_dB, ...
    'MajorTicks',-2:2:20,'ValueChangedFcn',@(s,~) setSNR(s.Value));
sSNR.Layout.Column = [1 2];

uilabel(gL,'Text','Jammer amp');
sJamA = uislider(gL,'Limits',[0 5],'Value',state.jammer_amp, ...
    'MajorTicks',0:0.5:5,'ValueChangedFcn',@(s,~) setJamAmp(s.Value));
sJamA.Layout.Column = [1 2];

uilabel(gL,'Text','Jammer Bandı Seç (B1..B5)'); uilabel(gL,'Text',' ');
jamBtnPanel = uipanel(gL,'BorderType','none'); jamBtnPanel.Layout.Column = [1 2];
gJB = uigridlayout(jamBtnPanel,[1,5],'ColumnWidth',{'1x','1x','1x','1x','1x'},'RowHeight',{'fit'});
btnJ = gobjects(1,5);
for k=1:5
    btnJ(k) = uibutton(gJB,'Text',sprintf('B%d (%.3f GHz)',k,state.bandsRF(k)/1e9), ...
        'ButtonPushedFcn',@(b,~) selectJammer(k));
end

cbWander = uicheckbox(gL,'Text','Jammer otomatik band gezer','Value',state.jamAutoWander, ...
    'ValueChangedFcn',@(c,~) setWander(c.Value));
cbWander.Layout.Column = [1 2];

uilabel(gL,'Text','Chunk size (bytes)');
eChunk = uieditfield(gL,'numeric','Limits',[1 Inf],'Value',state.chunk_size, ...
    'RoundFractionalValues','on','ValueChangedFcn',@(e,~) setChunkSize(e.Value));
eChunk.Layout.Column = [1 2];

btnStart = uibutton(gL,'Text','BAŞLAT','BackgroundColor',[0.21 0.6 0.24], ...
    'FontWeight','bold','ButtonPushedFcn',@(b,~) onStartStop());
btnStop  = uibutton(gL,'Text','DURDUR','BackgroundColor',[0.75 0.2 0.2], ...
    'Enable','off','ButtonPushedFcn',@(b,~) onStop());
btnStop.Layout.Column = 2;

lblVideo = uilabel(gL,'Text','Seçili video: yok','WordWrap','on'); lblVideo.Layout.Column = [1 2];
lblStats = uilabel(gL,'Text','Durum: hazır','WordWrap','on');      lblStats.Layout.Column = [1 2];

pbOuter = uipanel(gL,'Title','İlerleme','BackgroundColor',[0.95 0.95 0.95]); pbOuter.Layout.Column = [1 2];
pbInner = uipanel(pbOuter,'BackgroundColor',[0.2 0.6 0.9],'Position',[3 3 0 16]);
lblProg = uilabel(gL,'Text','0%'); lblProg.Layout.Column = [1 2];

% Mid plots
pMid = uipanel(gl,'Title','Grafikler','FontWeight','bold'); pMid.Layout.Row = [1 2]; pMid.Layout.Column = 2;
gM = uigridlayout(pMid,[3,1],'RowHeight',{'1x','1x','1x'},'Padding',6,'RowSpacing',8);
axSpec  = uiaxes(gM); title(axSpec,'Spektrum (IF, tıkla: en yakın banda snap)'); xlabel(axSpec,'Frekans (IF, Hz)'); ylabel(axSpec,'PSD');
axSpec.HitTest = 'on'; axSpec.PickableParts = 'all'; axSpec.ButtonDownFcn = @(~,~) onSpectrumClick();
axConst = uiaxes(gM); title(axConst,'Konstelasyon'); xlabel(axConst,'I'); ylabel(axConst,'Q');
axLLR   = uiaxes(gM); title(axLLR,'LLR Histogramı'); xlabel(axLLR,'LLR'); ylabel(axLLR,'Adet');

% Right metrics
pRight = uipanel(gl,'Title','Metrikler & Kayıt','FontWeight','bold');
pRight.Layout.Row = [1 2]; pRight.Layout.Column = 3;
gR = uigridlayout(pRight,[6,1],'RowHeight',{'fit','fit','1x','fit','1x','1x'},'Padding',8,'RowSpacing',10);

gTop = uigridlayout(gR,[3,2],'RowHeight',{'fit','fit','fit'},'ColumnWidth',{'1x','1x'},'RowSpacing',6);
uilabel(gTop,'Text','SNR (dB):');         lblSNR = uilabel(gTop,'Text',num2str(state.EbNo_dB,'%0.2f'));
uilabel(gTop,'Text','Çalışma CF (GHz):'); lblCF  = uilabel(gTop,'Text',num2str(state.bandsRF(state.cfIdx)/1e9,'%.3f'));
uilabel(gTop,'Text','Jammer f (GHz):');   lblJam = uilabel(gTop,'Text',num2str(state.bandsRF(state.jamIdx)/1e9,'%.3f'));

gBER = uigridlayout(gR,[2,2],'RowHeight',{'fit','fit'},'ColumnWidth',{'1x','1x'});
uilabel(gBER,'Text','BER (öncesi):');  lblBERpre  = uilabel(gBER,'Text','-');
uilabel(gBER,'Text','BER (sonrası):'); lblBERpost = uilabel(gBER,'Text','-');

axBER = uiaxes(gR); title(axBER,'BER Zaman Serisi'); xlabel(axBER,'Parça #'); ylabel(axBER,'BER');
hold(axBER,'on'); pltPre = plot(axBER,nan,nan,'-'); pltPost = plot(axBER,nan,nan,'-'); legend(axBER,{'Pre-LDPC','Post-LDPC'},'Location','northeast');

uilabel(gR,'Text','5 Bant Taraması');
tbl = uitable(gR,'ColumnName',{'Band','RF_GHz','Power_dB','BG_dB','Thresh','Durum'},'Data',cell(0,6));

ta = uitextarea(gR,'Editable','off','FontName','Consolas');

lblHelp = uilabel(gl,'Text','İpucu: Jammer butonları RF (GHz). Spektrumda tıklayınca IF’te en yakın banda snap eder; metrikler RF gösterir.');
lblHelp.Layout.Row = 3; lblHelp.Layout.Column = [2 3];

%% -------------------- CALLBACKS --------------------
    function onChooseVideo()
        [f,p] = uigetfile({'*.mp4;*.bin;*.dat','Video/Raw (*.mp4,*.bin,*.dat)';'*.*','Hepsi'});
        if isequal(f,0), return; end
        state.videoPath = fullfile(p,f);
        lblVideo.Text = "Seçili video: " + f;
        logmsg("Video seçildi: " + state.videoPath);
    end

    function onChooseOutput()
        p = uigetdir();
        if isequal(p,0), return; end
        state.outputFolder = p;
        lblOut.Text = "Çıkış: " + p;
        logmsg("Çıkış klasörü: " + p);
    end

    function setSNR(v)
        state.EbNo_dB = v; lblSNR.Text = sprintf('%.2f',v);
    end

    function setJamAmp(v)
        state.jammer_amp = v;
    end

    function setChunkSize(v)
        v = max(1000,round(v)); state.chunk_size = v; eChunk.Value = v;
    end

    function setWander(tf)
        state.jamAutoWander = tf;
        % Oto gez açıkken butonları pasifleştir
        for k=1:5
            if tf, btnJ(k).Enable = 'off'; else, btnJ(k).Enable = 'on'; end
        end
    end

    function selectJammer(k)
        state.jamIdx = k;
        lblJam.Text = sprintf('%.3f', state.bandsRF(k)/1e9);
        logmsg(sprintf("Jammer band: B%d (%.3f GHz)",k,state.bandsRF(k)/1e9));
    end

    function onSpectrumClick()
        cp = axSpec.CurrentPoint; x = cp(1,1); % IF Hz
        [~,idx] = min(abs(state.bandsIF - x));
        state.cfIdx = idx; lblCF.Text = sprintf('%.3f',state.bandsRF(idx)/1e9);
        logmsg(sprintf("CF grafikten: B%d (%.3f GHz)",idx,state.bandsRF(idx)/1e9));
    end

    function onStartStop()
        if state.running, return; end
        if state.videoPath == "", uialert(fig,'Lütfen bir video dosyası seçin.','Uyarı'); return; end
        state.running = true; state.stopRequested = false;
        set(btnStart,'Enable','off'); set(btnStop,'Enable','on');
        lblStats.Text = 'Durum: çalışıyor...';
        runSimulation();
        state.running = false; set(btnStart,'Enable','on'); set(btnStop,'Enable','off');
        lblStats.Text = 'Durum: durdu / bitti';
    end

    function onStop()
        state.stopRequested = true; lblStats.Text = 'Durum: durdurma isteniyor...';
    end

%% -------------------- CORE --------------------
    function runSimulation()
        [video_bytes,total_bytes] = readBytes(state.videoPath);
        chunk_size  = state.chunk_size;
        num_chunks  = ceil(double(total_bytes)/double(chunk_size));
        state.all_rx_bytes = uint8([]);
        berPreSeries  = nan(1,num_chunks); berPostSeries = nan(1,num_chunks);

        outFile = resolveOutFile();

        for part = 1:num_chunks
            if state.stopRequested, logmsg("Durdurma alındı."); break; end

            EbNo_dB     = state.EbNo_dB;
            jammer_amp  = state.jammer_amp;

            % Jammer bandı (IF)
            if state.jamAutoWander
                state.jamIdx = randi([1 5]);
                lblJam.Text = sprintf('%.3f', state.bandsRF(state.jamIdx)/1e9);
            end
            jammer_freq_if = state.bandsIF(state.jamIdx);

            % Çalışma bandı (IF) - mevcut seçim
            cf_if = state.bandsIF(state.cfIdx);

            % Parçayı al
            idx_start = (part-1)*chunk_size + 1;
            idx_end   = min(part*chunk_size, total_bytes);
            part_bytes = video_bytes(idx_start:idx_end);

            % Byte -> bit
            video_bits = reshape(de2bi(part_bytes,8,'left-msb')',[],1);

            % CRC ekle
            txCRC = state.crcGen(logical(video_bits));
            % Padding
            pad = mod(state.K - mod(length(txCRC),state.K), state.K);
            if pad>0, txCRC = [txCRC; false(pad,1)]; end

            % LDPC encode
            numBl = length(txCRC)/state.K;
            encBits = false(numBl*state.N,1);
            for i=1:numBl
                encBits((i-1)*state.N+1:i*state.N) = ldpcEncode(txCRC((i-1)*state.K+1:i*state.K),encCfg);
            end

            % 5 bant taraması ve kaçış (IF üzerinde ölç; tabloda RF yaz)
            [clean_if, tableRows] = scanJammer5(jammer_freq_if, jammer_amp, state.bandsIF, state.bandsRF, state.fs, state.T, EbNo_dB);
            tbl.Data = tableRows;
            if ~isempty(clean_if)
                [~,idx] = min(abs(state.bandsIF - clean_if(1)));
                state.cfIdx = idx; cf_if = state.bandsIF(idx);
                lblCF.Text = sprintf('%.3f', state.bandsRF(idx)/1e9);
            end

            % Görsel kanal (IF)
            t = 0:1/state.fs:state.T-1/state.fs;
            signal = sin(2*pi*cf_if*t);
            jammer = jammer_amp*sin(2*pi*jammer_freq_if*t);
            N0_vis = 1/(10^(EbNo_dB/10));
            noise_vis = sqrt(N0_vis/2)*randn(size(t));
            rx_total = signal + jammer + noise_vis;

            % Spektrum (IF ekseni, RF etiketli çizgiler)
            plotPSD5(axSpec, rx_total, state.fs, state.bandsIF, state.bandsRF, cf_if, jammer_freq_if);

            % Modülasyon & demodülasyon (Eb/N0'a göre)
            [llr, rxSym, symbols, modType] = modChanDemod(encBits, EbNo_dB);
            plotConst(axConst, symbols, rxSym, modType, EbNo_dB);
            plotLLR(axLLR, llr);

            % BER
            rxHard = llr(:) < 0;
            ber_pre = mean(rxHard ~= encBits);   lblBERpre.Text  = sprintf('%.6f', ber_pre);

            decBits = false(numBl*state.K,1);
            for i=1:numBl
                decBits((i-1)*state.K+1:i*state.K) = ldpcDecode(llr((i-1)*state.N+1:i*state.N), decCfg, state.maxItersLDPC);
            end
            [rxCRC, err] = state.crcDet(decBits);
            ber_post = mean(decBits(1:length(txCRC)) ~= txCRC); lblBERpost.Text = sprintf('%.6f', ber_post);
            if err
                logmsg(sprintf("Parça %d: CRC başarısız!", part));
                rxCRC = decBits(1:end - pad);
            else
                if pad>0 && length(rxCRC)>=pad, rxCRC = rxCRC(1:end-pad); end
            end

            % Bit -> byte ve biriktir
            padbits = mod(8 - mod(length(rxCRC),8), 8);
            if padbits>0, rxCRC = [rxCRC; false(padbits,1)]; end
            video_bytes_rx = uint8(bi2de(reshape(rxCRC,8,[])','left-msb'));
            state.all_rx_bytes = [state.all_rx_bytes; video_bytes_rx]; %#ok<AGROW>

            % Seriler & ilerleme
            berPreSeries(part)  = ber_pre;  berPostSeries(part) = ber_post;
            set(pltPre,'XData',1:part,'YData',berPreSeries(1:part));
            set(pltPost,'XData',1:part,'YData',berPostSeries(1:part));
            state.progress = double(idx_end)/double(total_bytes);
            updateProgress(state.progress);

            lblStats.Text = sprintf('Parça %d/%d | Mod: %s | CF=%.3f GHz | Jam=%.3f GHz', ...
                                    part,num_chunks,modType, ...
                                    state.bandsRF(state.cfIdx)/1e9, state.bandsRF(state.jamIdx)/1e9);
            drawnow limitrate;
        end

        % Çıktı
        if ~isempty(state.all_rx_bytes)
            try
                fid = fopen(outFile,'wb'); fwrite(fid,state.all_rx_bytes,'uint8'); fclose(fid);
                logmsg("Çıktı yazıldı: " + outFile);
            catch ME
                logmsg("Çıkış yazım hatası: " + ME.message);
            end
        end
    end

%% -------------------- HELPERS --------------------
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
            logmsg(sprintf("Video boyutu: %d byte", total_bytes));
        catch ME %#ok<NASGU>
            logmsg("Dosya okuma hatası. Sentetik test verisi oluşturuluyor (1 MB)...");
            video_bytes = uint8(randi([0 255],1e6,1)); total_bytes = numel(video_bytes);
        end
    end

    function [clean_centers_if, tableRows] = scanJammer5(jam_f_if, jam_amp, bandsIF, bandsRF, fs, T, EbNo_dB)
        % IF üzerinde jammer gücünü her band civarında ölç, RF etiketle
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

            margin = 500; % Hz
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
        PSD = abs(fft(rx_total,Nfft)/Nfft).^2;
        PSD = PSD(1:Nfft/2+1);

        plot(ax, f, 10*log10(PSD+eps), 'LineWidth', 1); grid(ax,'on'); hold(ax,'on');
        yl = ylim(ax);
        for k=1:numel(bandsIF)
            xline(ax, bandsIF(k), '-', sprintf('B%d (%.3f GHz)', k, bandsRF(k)/1e9), ...
                'LineWidth',1.0,'Alpha',0.6);
        end
        xline(ax, cf_if, 'g-',  'LineWidth',1.6);
        xline(ax, jf_if, 'r--', 'LineWidth',1.6);
        ylim(ax,yl);
        legend(ax,{'PSD','Bands (RF labels)','CF(IF)','Jam(IF)'},'Location','northeastoutside');
        set(findobj(ax,'Type','ConstantLine'),'HitTest','off','PickableParts','none');
        hold(ax,'off');
        xlabel(ax,'Frekans (IF, Hz)'); ylabel(ax,'PSD');
        title(ax,'Spektrum (IF) – RF etiketli bantlar');
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
        legend(ax,{'Tx','Rx'},'Location','northeastoutside');
        hold(ax,'off');
    end

    function plotLLR(ax, llr)
        cla(ax); histogram(ax, llr, 100); grid(ax,'on');
        xlabel(ax,'LLR'); ylabel(ax,'Adet'); title(ax,'LLR Dağılımı');
    end

    function updateProgress(p)
        p = max(0,min(1,p));
        lblProg.Text = sprintf('%.1f%%', 100*p);
        outerW = pbOuter.Position(3);
        set(pbInner,'Position',[3 3 max(0,(outerW-6)*p) 16]);
    end

    function logmsg(s)
        ta.Value = [ta.Value; string(datestr(now,'HH:MM:SS')) + " | " + string(s)];
        drawnow limitrate;
    end

    function m = select_modulation(EbNo_dB)
        if EbNo_dB < 4
            m = 'BPSK';
        elseif EbNo_dB < 8
            m = '4QAM';
        elseif EbNo_dB < 15
            m = '16QAM';
        else
            m = '64QAM';
        end
    end
end
