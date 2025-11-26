function pluto_tx_1
% LTE-like 1.4 MHz OFDM TX (PSS+SSS+CRS+Data), Pluto + GUI (R2020+ uyumlu)
% Grid: Nfft=128, Fs=1.92e6 (Δf=15kHz), normal CP [10,9,9,9,9,9,9]
% NOT: Gerçek LTE ile interoperable değildir ama RX ile senkron çalışır.

    app = buildUI();
    app.isRunning = false;
    guidata(app.fig, app);

function app = buildUI()
    app.fig = uifigure('Name','Pluto TX (LTE-like 1.4 MHz)','Position',[100 100 860 560]);
    g = uigridlayout(app.fig,[6 6]);
    g.RowHeight  = {28,28,28,28,28,'1x'};
    g.ColumnWidth= {140,140,140,140,120,'1x'};

    % --- TX Serial ---
    lbl = uilabel(g,'Text','TX Serial','HorizontalAlignment','right');
    lbl.Layout.Row = 1; lbl.Layout.Column = 1;
    app.edSerial = uieditfield(g,'text','Value','sn:1044735411960005efff21000d2d1fa1f5');
    app.edSerial.Layout.Row = 1; app.edSerial.Layout.Column = [2 4];

    % --- Fc ---
    lbl = uilabel(g,'Text','Fc [Hz]','HorizontalAlignment','right');
    lbl.Layout.Row = 1; lbl.Layout.Column = 5;
    app.edFc = uieditfield(g,'numeric','Value',915e6,'Limits',[70e6 6e9],'ValueDisplayFormat','%.0f');
    app.edFc.Layout.Row = 1; app.edFc.Layout.Column = 6;

    % --- Gain ---
    lbl = uilabel(g,'Text','TX Gain [dB]','HorizontalAlignment','right');
    lbl.Layout.Row = 2; lbl.Layout.Column = 1;
    app.sldGain = uislider(g,'Limits',[-60 0],'Value',-10);
    app.sldGain.Layout.Row = 2; app.sldGain.Layout.Column = [2 6];

    % --- Fs ---
    lbl = uilabel(g,'Text','Fs','HorizontalAlignment','right');
    lbl.Layout.Row = 3; lbl.Layout.Column = 1;
    app.ddFs = uidropdown(g,'Items',{'1.92e6','3.84e6'},'Value','1.92e6');
    app.ddFs.Layout.Row = 3; app.ddFs.Layout.Column = 2;

    % --- CellID ---
    lbl = uilabel(g,'Text','CellID (0..503)','HorizontalAlignment','right');
    lbl.Layout.Row = 3; lbl.Layout.Column = 3;
    app.edCellID = uieditfield(g,'numeric','Value',0,'Limits',[0 503]);
    app.edCellID.Layout.Row = 3; app.edCellID.Layout.Column = 4;

    % --- XOR key-stream (opsiyonel) ---
    app.chkCipher = uicheckbox(g,'Text','XOR key-stream (opsiyonel)','Value',false);
    app.chkCipher.Layout.Row = 3; app.chkCipher.Layout.Column = [5 6];

    % --- Payload seç ---
    app.btnFile = uibutton(g,'Text','Payload Dosyası Seç','ButtonPushedFcn',@pickFile);
    app.btnFile.Layout.Row = 4; app.btnFile.Layout.Column = 1;
    app.labFile = uilabel(g,'Text','[opsiyonel]');
    app.labFile.Layout.Row = 4; app.labFile.Layout.Column = [2 6];

    % --- Start/Stop ---
    app.btnStart = uibutton(g,'Text','BAŞLAT','ButtonPushedFcn',@startTx,'BackgroundColor',[0.2 0.7 0.2]);
    app.btnStart.Layout.Row = 5; app.btnStart.Layout.Column = 1;
    app.btnStop  = uibutton(g,'Text','DURDUR','ButtonPushedFcn',@stopTx,'BackgroundColor',[0.8 0.2 0.2]);
    app.btnStop.Layout.Row = 5; app.btnStop.Layout.Column = 2;

    % --- Spektrum ---
    app.axSpec = uiaxes(g);
    app.axSpec.Layout.Row = 6; app.axSpec.Layout.Column = [1 6];
    title(app.axSpec,'TX Spectrum (son sembol)'); xlabel(app.axSpec,'FFT bin'); ylabel(app.axSpec,'|X|');

    % State
    app.payloadPath = '';
    app.txObj = [];
end


    function pickFile(~,~)
        [f,p] = uigetfile({'*.*','All Files'},'Payload dosyası');
        if isequal(f,0), return; end
        app = guidata(app.fig); app.payloadPath = fullfile(p,f);
        app.labFile.Text = app.payloadPath; guidata(app.fig,app);
    end

    function startTx(~,~)
        app = guidata(app.fig);
        if app.isRunning, return; end
        app.isRunning = true; guidata(app.fig,app);

        % --- Parametreler ---
        Fc   = app.edFc.Value;
        Fs   = str2double(app.ddFs.Value);
        gain = app.sldGain.Value;
        radioID = strtrim(app.edSerial.Value);
        NcellID = app.edCellID.Value;
        useCipher = app.chkCipher.Value;

        % Numerology
        Nfft = 128;
        scs  = Fs/Nfft; % 15 kHz hedefi
        if abs(scs-15e3) > 1e-6
            uialert(app.fig,'Fs/Nfft = 15kHz değil; Fs=1.92e6 seçin (Nfft=128).','Uyarı');
        end
        cpSym = [10 9 9 9 9 9 9]; % normal CP
        RB = 6; Nsc = RB*12;     % 72 aktif taşıyıcı
        dc = Nfft/2 + 1;
        kL = dc-(Nsc/2):dc-1; kR = dc+1:dc+(Nsc/2);
        actIdx = [kL kR];

        % Pluto init
        try
            Tx = sdrtx('Pluto','RadioID',radioID,'CenterFrequency',Fc,'BasebandSampleRate',Fs,'Gain',gain);
        catch ME
            uialert(app.fig,ME.message,'Pluto Hatası'); stopTx(); return;
        end
        app.txObj = Tx; guidata(app.fig,app);

        % Payload
        if ~isempty(app.payloadPath)
            fid = fopen(app.payloadPath,'rb'); 
            if fid<0, uialert(app.fig,'Dosya açılamadı.','Hata'); stopTx(); return; end
            fileBytes = fread(fid,'*uint8'); fclose(fid);
        else
            fileBytes = [];
        end
        filePtr = 1;

        rngScr = @(sf,sz) gold_seq( (NcellID*2^9 + sf), sz );
        keyBytes = uint8([]);
        if useCipher, keyBytes = uint8(prbs_bytes(4096)); end

        sf = 0; % subframe 0..9
        try
            while true
                app = guidata(app.fig);
                if ~app.isRunning, break; end

                [timeBlock, Xlast] = buildSubframe(sf, Nfft, cpSym, actIdx, NcellID, ...
                    @pullBits, rngScr, keyBytes);

                Tx(timeBlock);

                if isvalid(app.axSpec)
                    stem(app.axSpec,abs(Xlast)); xlim(app.axSpec,[1 Nfft]); drawnow limitrate;
                end

                sf = mod(sf+1,10);
            end
        catch ME
            uialert(app.fig,ME.message,'TX Döngü Hatası');
        end

        release(Tx); app.txObj = []; guidata(app.fig,app);

        function bits = pullBits(nbits)
            if isempty(fileBytes)
                bits = prbs_bits(nbits);
            else
                nBytesNeed = ceil(nbits/8);
                if filePtr + nBytesNeed - 1 > numel(fileBytes)
                    left = numel(fileBytes)-filePtr+1;
                    pool = [fileBytes(filePtr:end); fileBytes(1:(nBytesNeed-left))];
                    filePtr = nBytesNeed-left+1;
                else
                    pool = fileBytes(filePtr:filePtr+nBytesNeed-1);
                    filePtr = filePtr + nBytesNeed;
                end
                bits = de2bi(pool,8,'left-msb').'; bits = bits(:);
                bits = bits(1:nbits);
            end
        end
    end

    function stopTx(~,~)
        app = guidata(app.fig);
        app.isRunning = false;
        if ~isempty(app.txObj)
            try, release(app.txObj); catch, end
            app.txObj = [];
        end
        guidata(app.fig,app);
    end
end

% --------- Subframe Builder (TX) ---------
function [y_time, Xlast] = buildSubframe(sf, Nfft, cpSym, actIdx, NcellID, pullBits, rngScr, keyBytes)
    NsymSF = 14;
    Nsc = numel(actIdx);
    M = 4; bps = log2(M);

    X = zeros(Nfft, NsymSF);
    dmMask = false(Nfft, NsymSF);

    % CRS-like pilots: semboller {1,5,8,12} (1-based)
    pilotSyms = [1,5,8,12];
    pilotPos  = actIdx(1:6:end);
    for s = pilotSyms
        X(pilotPos, s) = crs_seq(length(pilotPos), NcellID, s-1);
        dmMask(pilotPos, s) = true;
    end

    % SSS/PSS subframe 0&5: slot1 sym5=12, sym6=13 (toplam 14 sembol var)
    if (sf==0 || sf==5)
        [kL,kR] = center62(actIdx);
        sss = sss_lite(NcellID);              % 62
        X(kL,12) = sss(1:31); X(kR,12) = sss(32:62); dmMask(kL,12)=true; dmMask(kR,12)=true;
        roots = [25 29 34]; r = roots(mod(NcellID,3)+1);
        pss = pss_zc_62(r);
        X(kL,13) = pss(1:31); X(kR,13) = pss(32:62); dmMask(kL,13)=true; dmMask(kR,13)=true;
    end

    % Data RE mask
    dataMask = false(Nfft, NsymSF); dataMask(actIdx,:)=true; dataMask(dmMask)=false;

    % Bit üret (scramble + opsiyonel XOR)
    nDataRE = nnz(dataMask); nbits = nDataRE*bps;
    bits = pullBits(nbits);
    if ~isempty(keyBytes)
        kbits = de2bi(repmat(keyBytes, ceil(nbits/8),1),8,'left-msb').'; kbits = kbits(:); kbits = kbits(1:nbits);
        bits = xor(bits, kbits);
    end
    c = rngScr(sf, nbits); bits = xor(bits, c);

    % QPSK
    symIdx = bi2de(reshape(bits,bps,[]).','left-msb');
    dQ = pskmod(symIdx,4,0,'gray');
    Xi = find(dataMask); X(Xi) = dQ;

    % OFDM + CP
    cpVec = repmat(cpSym,1,2); y = [];
    for s=1:NsymSF
        xt = ifft(ifftshift(X(:,s)));
        cpN = cpVec(s);
        xt = [xt(end-cpN+1:end); xt];
        y  = [y; xt]; %#ok<AGROW>
    end
    y_time = single(y ./ max(1e-6, max(abs(y))));
    Xlast = abs(X(:,end));
end

% --------- Ortak yardımcılar ---------
function [kL,kR] = center62(actIdx)
    mid = actIdx(ceil(numel(actIdx)/2));
    leftAll = actIdx(actIdx<mid);
    rightAll= actIdx(actIdx>mid);
    kL = leftAll(end-31+1:end);
    kR = rightAll(1:31);
end

function s = crs_seq(N, cellID, sym)
    b = gold_seq(cellID*1024 + sym, 2*N);
    b = reshape(b,2,[]); idx = bi2de(b.','left-msb');
    s = pskmod(idx,4,0,'gray');
end

function pss = pss_zc_62(u)
    Nzc = 63; n=(0:Nzc-1).';
    zc = exp(-1j*pi*u.*n.*(n+1)/Nzc);
    zc((Nzc+1)/2) = []; 
    pss = zc / sqrt(mean(abs(zc).^2));
end

function sss = sss_lite(cellID)
    % Basitleştirilmiş SSS (62 uzunluk) – demoda yeterli
    N = 62;
    x1 = mseq_poly(7,[7 1],1,N);
    x2 = mseq_poly(7,[7 3],mod(cellID,127)+1,N);
    b = xor(x1,x2);
    sss = 2*double(b)-1; sss = sss(:);
end

function bits = mseq_poly(m,taps,ini,N)
    reg = zeros(1,m); reg(ini)=1;
    out = zeros(1,N);
    for i=1:N
        out(i) = reg(end);
        fb = 0;
        for t=taps, fb = xor(fb, reg(m - t + 1)); end
        reg = [fb reg(1:end-1)];
    end
    bits = out;
end

function b = gold_seq(cinit,N)
    b = zeros(N,1);
    x1 = ones(31,1);
    x2 = de2bi(mod(cinit,2^31),31,'left-msb').'; x2(x2==0)=1;
    for n=1:N
        b(n) = xor(x1(end), x2(end));
        x1fb = xor(x1(1), x1(4));        % (31,28) eşleniği basitçe
        x1 = [x1fb; x1(1:end-1)];
        x2fb = xor(x2(4), xor(x2(5), xor(x2(6), x2(21))));
        x2 = [x2fb; x2(1:end-1)];
    end
    b = logical(b);
end

function b = prbs_bits(N)
    s = uint32(0xACE1);
    b = false(N,1);
    for i=1:N
        s = bitxor(bitshift(s,1,'uint32'), bitand(bitxor(s, bitshift(s,-2)),1));
        b(i) = bitand(s,1);
    end
end

function by = prbs_bytes(N)
    b = prbs_bits(8*N);
    bb = reshape(b,8,[]).';
    by = uint8(bi2de(bb,'left-msb'));
end
