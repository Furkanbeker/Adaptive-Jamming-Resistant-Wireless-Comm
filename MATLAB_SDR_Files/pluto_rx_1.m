function pluto_rx_1
% LTE-like 1.4 MHz OFDM RX (Pluto + GUI)
% PSS tabanlı kaba zaman/CFO, CP hizalama, FFT, CRS eşitleme, QPSK demod.

    app = buildUI();  app.isRunning=false; guidata(app.fig,app);

function app = buildUI()
    app.fig = uifigure('Name','Pluto RX (LTE-like 1.4 MHz)','Position',[100 100 900 640]);
    g = uigridlayout(app.fig,[7 6]);
    g.RowHeight  = {28,28,28,28,28,220,'1x'};
    g.ColumnWidth= {140,140,140,140,120,'1x'};

    % --- RX Serial ---
    lbl = uilabel(g,'Text','RX Serial','HorizontalAlignment','right');
    lbl.Layout.Row = 1; lbl.Layout.Column = 1;
    app.edSerial = uieditfield(g,'text','Value','sn:10447354119600060b003a0007241ed971');
    app.edSerial.Layout.Row = 1; app.edSerial.Layout.Column = [2 4];

    % --- Fc ---
    lbl = uilabel(g,'Text','Fc [Hz]','HorizontalAlignment','right');
    lbl.Layout.Row = 1; lbl.Layout.Column = 5;
    app.edFc = uieditfield(g,'numeric','Value',915e6,'Limits',[70e6 6e9],'ValueDisplayFormat','%.0f');
    app.edFc.Layout.Row = 1; app.edFc.Layout.Column = 6;

    % --- Gain ---
    lbl = uilabel(g,'Text','RX Gain [dB]','HorizontalAlignment','right');
    lbl.Layout.Row = 2; lbl.Layout.Column = 1;
    app.sldGain = uislider(g,'Limits',[0 60],'Value',30);
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

    % --- Start/Stop ---
    app.btnStart = uibutton(g,'Text','BAŞLAT','ButtonPushedFcn',@startRx,'BackgroundColor',[0.2 0.7 0.2]);
    app.btnStart.Layout.Row = 4; app.btnStart.Layout.Column = 1;
    app.btnStop  = uibutton(g,'Text','DURDUR','ButtonPushedFcn',@stopRx,'BackgroundColor',[0.8 0.2 0.2]);
    app.btnStop.Layout.Row = 4; app.btnStop.Layout.Column = 2;

    % --- Info ---
    app.labInfo = uilabel(g,'Text','Hazır.');
    app.labInfo.Layout.Row = 4; app.labInfo.Layout.Column = [3 6];

    % --- Plots ---
    app.axConst = uiaxes(g);
    app.axConst.Layout.Row = 6; app.axConst.Layout.Column = [1 3];
    title(app.axConst,'Constellation'); xlabel(app.axConst,'I'); ylabel(app.axConst,'Q');

    app.axSpec  = uiaxes(g);
    app.axSpec.Layout.Row = 6; app.axSpec.Layout.Column = [4 6];
    title(app.axSpec,'PSS Corr / Spectrum');

    % --- Log ---
    app.txtLog = uitextarea(g,'Editable','off');
    app.txtLog.Layout.Row = 7; app.txtLog.Layout.Column = [1 6];

    app.rxObj = [];
end


    function startRx(~,~)
        app = guidata(app.fig);
        if app.isRunning, return; end
        app.isRunning = true; guidata(app.fig,app);

        Fc   = app.edFc.Value;
        Fs   = str2double(app.ddFs.Value);
        gain = app.sldGain.Value;
        radioID = strtrim(app.edSerial.Value);
        NcellID = app.edCellID.Value;

        % Numerology
        Nfft = 128;
        Ts = 1/Fs;
        cpSym = [10 9 9 9 9 9 9];
        RB=6; Nsc=RB*12;
        dc = Nfft/2 + 1;
        kL = dc-(Nsc/2):dc-1; kR = dc+1:dc+(Nsc/2);
        actIdx = [kL kR];

        % PSS set
        uVals = [25 29 34];
        PSS = cell(1,3);
        for i=1:3, PSS{i} = pss_zc_62(uVals(i)); end
        [k62L,k62R] = center62(actIdx);

        % Pluto init
        spf = 4096;
        try
            Rx = sdrrx('Pluto','RadioID',radioID,'CenterFrequency',Fc,'BasebandSampleRate',Fs,...
                'OutputDataType','double','GainSource','Manual','Gain',gain,'SamplesPerFrame',spf,'EnableBurstMode',false);
        catch ME
            uialert(app.fig,ME.message,'Pluto Hatası'); stopRx(); return;
        end
        app.rxObj = Rx; guidata(app.fig,app);

        cfoEst = 0; lock = false; evmEMA = 0.0;
        while true
            app = guidata(app.fig);
            if ~app.isRunning, break; end

            y = Rx(); if isempty(y), pause(0.001); continue; end

            % CFO düzelt (mevcut kestirimle)
            if cfoEst~=0
                n = (0:numel(y)-1).';
                y = y .* exp(-1j*2*pi*cfoEst*n/Fs);
            end

            % PSS araması (basit matched)
            corrVal=0; bestIdx=1; bestU=1;
            if numel(y) >= (Nfft+10)
                for i=1:3
                    X0 = zeros(Nfft,1); X0(k62L)=PSS{i}(1:31); X0(k62R)=PSS{i}(32:62);
                    s0 = ifft(ifftshift(X0));
                    c = abs(conv(y, flipud(conj(s0)), 'valid'));
                    [mv,mi]= max(c);
                    if mv>corrVal, corrVal=mv; bestIdx=mi; bestU=i; end
                end
            end

            % Eşik ve CFO
            if corrVal>50
                lock = true;
                L = Nfft;
                if bestIdx+L*2 <= numel(y)
                    seg1 = y(bestIdx:bestIdx+L-1);
                    seg2 = y(bestIdx+L:bestIdx+2*L-1);
                    ang = angle(sum(conj(seg1).*seg2));
                    cfoEst = ang/(2*pi*L*Ts);
                end
            end

            if lock
                [y_sf, ok] = grab_subframe(y, bestIdx, Nfft, cpSym);
                if ~ok, continue; end

                [Xsf, ~] = ofdm_fft(y_sf, Nfft, cpSym);

                pilotSyms = [1,5,8,12];
                pilotPos  = actIdx(1:6:end);
                Hhat = ones(numel(actIdx),14);
                for s = pilotSyms
                    ps = crs_seq(length(pilotPos), NcellID, s-1);
                    Yp = Xsf(pilotPos, s);
                    Hs = Yp ./ ps;
                    Hhat(:,s) = interp1(pilotPos.', Hs, actIdx.','linear','extrap');
                end
                H = mean(Hhat(:,pilotSyms),2); H = repmat(H,1,14);

                Xeq = Xsf;
                for s=1:14
                    Xeq(actIdx,s) = Xsf(actIdx,s) ./ (H(:,s)+1e-8);
                end

                % Maskeler
                dmMask = false(Nfft,14);
                for s = pilotSyms, dmMask(pilotPos,s)=true; end
                dmMask(k62L,12)=true; dmMask(k62R,12)=true; % SSS
                dmMask(k62L,13)=true; dmMask(k62R,13)=true; % PSS
                dataMask = false(Nfft,14); dataMask(actIdx,:)=true; dataMask(dmMask)=false;

                Xi = find(dataMask);
                dQ = Xeq(Xi);

                % Constellation
                if isvalid(app.axConst)
                    plot(app.axConst, real(dQ), imag(dQ), '.'); axis(app.axConst,[-2 2 -2 2]); grid(app.axConst,'on');
                end

                % EVM
                dem = pskdemod(dQ,4,0,'gray');
                ref = pskmod(dem,4,0,'gray');
                evm = sqrt(mean(abs(dQ - ref).^2) / mean(abs(ref).^2));
                evmEMA = 0.9*evmEMA + 0.1*evm;

                % Log (BURASI DÜZELTİLDİ: literal diziye indeksleme yok)
                rootVal = uVals(bestU);
                msg = sprintf('Lock=1, PSSroot=%d, CFO=%.1f Hz, EVM=%.3f (EMA=%.3f), dRE=%d', ...
                               rootVal, cfoEst, evm, evmEMA, numel(dQ));
                app.txtLog.Value = [app.txtLog.Value; msg];
                app.txtLog.Value = app.txtLog.Value(max(1,end-120):end);

                if isvalid(app.axSpec)
                    Yplot = y(end-min(2048,numel(y))+1:end);
                    plot(app.axSpec, abs(fftshift(fft(Yplot))));
                    title(app.axSpec, sprintf('PSS corr=%.1f',corrVal));
                end
            end

            drawnow limitrate
        end

        release(Rx); app.rxObj=[]; guidata(app.fig,app);
    end

    function stopRx(~,~)
        app = guidata(app.fig);
        app.isRunning=false;
        if ~isempty(app.rxObj)
            try, release(app.rxObj); catch, end
            app.rxObj=[];
        end
        guidata(app.fig,app);
    end
end

% -------- Helpers --------
function [ysf, ok] = grab_subframe(y, idx, Nfft, cpSym)
    cpVec = repmat(cpSym,1,2);
    need = sum(Nfft+cpVec);
    if idx+need-1 > numel(y), ok=false; ysf=[]; return; end
    ysf = y(idx:idx+need-1); ok=true;
end

function [Xsf, ofs] = ofdm_fft(y_sf, Nfft, cpSym)
    cpVec = repmat(cpSym,1,2);
    Xsf = zeros(Nfft,14); ofs = zeros(1,14);
    p = 1;
    for s=1:14
        cpN = cpVec(s);
        seg = y_sf(p+cpN : p+cpN+Nfft-1);
        ofs(s) = p;
        Xsf(:,s) = fftshift(fft(seg));
        p = p + (Nfft+cpN);
    end
end

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
    Nzc=63; n=(0:Nzc-1).';
    zc = exp(-1j*pi*u.*n.*(n+1)/Nzc);
    zc((Nzc+1)/2) = []; pss = zc / sqrt(mean(abs(zc).^2));
end

function b = gold_seq(cinit,N)
    b = zeros(N,1);
    x1 = ones(31,1);
    x2 = de2bi(mod(cinit,2^31),31,'left-msb').'; x2(x2==0)=1;
    for n=1:N
        b(n) = xor(x1(end), x2(end));
        x1fb = xor(x1(1), x1(4));
        x1 = [x1fb; x1(1:end-1)];
        x2fb = xor(x2(4), xor(x2(5), xor(x2(6), x2(21))));
        x2 = [x2fb; x2(1:end-1)];
    end
    b = logical(b);
end
