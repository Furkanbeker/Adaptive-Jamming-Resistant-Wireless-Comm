clear; close all; clc;

channel = "OverTheAir";
deviceName = "Pluto";
centerFrequency = 1.8e9;
txGain = 0;

usbRx='sn:10447354119600060b003a0007241ed971';
usbTx='sn:1044735411960005efff21000d2d1fa1f5';

warning('off','plutoradio:sysobj:FirmwareIncompatible');

% Chose input filePath
% [filename, pathname] = uigetfile('*.txt', 'Select Input Text File');
% if isequal(filename,0)
%     error('File selection canceled.');
% end
% filePath = fullfile(pathname, filename);
% 
% fid = fopen(filePath, 'r');
% txt = fread(fid, '*char')';
% fclose(fid);
% 

filePath = 'C:\Users\mfb36\Desktop\input.txt';
fid = fopen(filePath, 'r');
txt = fread(fid, '*char')';
fclose(fid);

binData = dec2bin(txt, 8);         % Convert each character to 8-bit binary (ASCII)
trData = reshape((binData - '0').', [], 1);  % Generate serialized bitstream

txsim.RC = 'R.5';                   % RMC type (10 MHz bandwidth)
txsim.NCellID = 53;
txsim.NFrame = 101;
txsim.NumAntennas = 1;

rmc = lteRMCDL(txsim.RC);
trBlkSize = rmc.PDSCH.TrBlkSizes;
frame_factor = ceil(numel(trData) / sum(trBlkSize(:)));
frame_factor = frame_factor / 5 + 3;
txsim.TotFrames = ceil((numel(trData) / frame_factor) / sum(trBlkSize(:)));

rmc.NCellID = txsim.NCellID;
rmc.NFrame = txsim.NFrame;
rmc.TotSubframes = txsim.TotFrames * 10;
rmc.CellRefP = txsim.NumAntennas;
rmc.PDSCH.RVSeq = 0;
rmc.OCNGPDSCHEnable = "Off";
rmc.OCNGPDCCHEnable = "Off";

if rmc.CellRefP >= 2
    rmc.PDSCH.TxScheme = "TxDiversity";
    rmc.OCNGPDSCH.TxScheme = "TxDiversity";
else
    rmc.PDSCH.TxScheme = "Port0";
    rmc.OCNGPDSCH.TxScheme = "Port0";
end

rmc.PDSCH.NLayers = txsim.NumAntennas;

fprintf("\nGenerating LTE transmit waveform:\nPacked data into %d frame(s).\n\n", txsim.TotFrames);

eNodeBOutput = [];
for i = 1:frame_factor
    txsim.NFrame = 100 + i;
    rmc.NFrame = txsim.NFrame;

    startIdx = floor((i-1) * length(trData) / frame_factor) + 1;
    endIdx = floor(i * length(trData) / frame_factor);
    dataChunk = trData(startIdx:endIdx);

    [eNodeBOutput(:,i), ~, rmc] = lteRMCDLTool(rmc, dataChunk);
end

sdrTransmitter = sdrtx(deviceName, 'RadioID', usbTx);
sdrTransmitter.BasebandSampleRate = rmc.SamplingRate;
sdrTransmitter.CenterFrequency = centerFrequency;
sdrTransmitter.Gain = txGain;
sdrTransmitter.ChannelMapping = 1:txsim.NumAntennas;

sdrReceiver = sdrrx(deviceName, 'RadioID', usbRx);
sdrReceiver.BasebandSampleRate = rmc.SamplingRate;
sdrReceiver.CenterFrequency = centerFrequency;
sdrReceiver.OutputDataType = "double";
sdrReceiver.ChannelMapping = 1:txsim.NumAntennas;

framesPerCapture = txsim.TotFrames + 1;
captureTime = framesPerCapture * 10e-3;  % 10 ms per frame

power = zeros(5,1);
freqs = [centerFrequency - 6e6, ...
         centerFrequency - 3e6, ...
         centerFrequency, ...
         centerFrequency + 3e6, ...
         centerFrequency + 6e6];

for k = 1:5
    sdrReceiver.CenterFrequency = freqs(k);
    jammingCheck = capture(sdrReceiver, captureTime, "Seconds");
    power(k) = pow2db(bandpower(jammingCheck));
end

[~, index] = min(power);
centerFrequency = freqs(index);
sdrTransmitter.CenterFrequency = centerFrequency;
sdrReceiver.CenterFrequency = centerFrequency;

fprintf("\nSelected frequency for transmission: %.2f MHz\n", centerFrequency / 1e6);

fprintf("Transmitting frame %d at %.2f MHz...\n", 1, centerFrequency / 1e6);

% Normalize to prevent RF saturation
powerScaleFactor = 1; % power scale of signal
eNodeBOutput = eNodeBOutput .* (1 ./ max(abs(eNodeBOutput)) * powerScaleFactor);

% Transmit one frame
fullWaveform = eNodeBOutput(:);  % Tüm frame'leri sırayla tek vektörde birleştir
transmitRepeat(sdrTransmitter, fullWaveform);

% Receive signal
rxWaveform = sdrReceiver();

% Optional: Add AWGN noise (manual control)
addNoise = true;               
noisePower = 0;             

if addNoise
    noise = sqrt(noisePower/2) * (randn(size(rxWaveform)) + 1j * randn(size(rxWaveform)));
    rxWaveform = rxWaveform + noise;
    fprintf("Added manual AWGN noise with power: %.6f\n", noisePower);
end


% Frequency offset correction
frequencyOffset = lteFrequencyOffset(rmc, rxWaveform);
rxWaveform = lteFrequencyCorrect(rmc, rxWaveform, frequencyOffset);
fprintf("Frequency offset corrected by %.2f Hz\n", frequencyOffset);

% Cell search
cellSearch.SSSDetection = "PostFFT";
cellSearch.MaxCellCount = 1;
[NCellID, frameOffset] = lteCellSearch(rmc, rxWaveform, cellSearch);
rmc.NCellID = NCellID;
fprintf("Cell ID: %d | Frame offset: %d samples\n", NCellID, frameOffset);

% Remove timing offset
rxWaveform = rxWaveform(frameOffset+1:end);

% Truncate to full frames only
samplesPerFrame = 10e-3 * rmc.SamplingRate;
tailSamples = mod(length(rxWaveform), samplesPerFrame);
rxWaveform = rxWaveform(1:end-tailSamples);

% OFDM demodulation
rxGrid = lteOFDMDemodulate(rmc, rxWaveform);

% Channel estimation configuration
cec.PilotAverage = "UserDefined";
cec.FreqWindow = 9;
cec.TimeWindow = 9;
cec.InterpType = "Cubic";
cec.InterpWindow = "Centered";
cec.InterpWinSize = 3;

% Estimate channel
[hest, nest] = lteDLChannelEstimate(rmc, cec, rxGrid);

sfDims = lteResourceGridSize(rmc);
Lsf = sfDims(2);          
LFrame = 10 * Lsf;        
numFullFrames = floor(length(rxWaveform) / (10e-3 * rmc.SamplingRate));
samplesPerFrame = 10e-3 * rmc.SamplingRate;

rxDataFrame = zeros(sum(rmc.PDSCH.TrBlkSizes(:)), numFullFrames);
recFrames = zeros(numFullFrames,1);

for frame = 0:(numFullFrames - 1)
    rmc.NSubframe = 0;
    rxsf = rxGrid(:, frame * LFrame + (1:Lsf), :);
    hestsf = hest(:, frame * LFrame + (1:Lsf), :, :);
    
    % PBCH decode
    pbchIndices = ltePBCHIndices(rmc);
    [pbchRx, pbchHest] = lteExtractResources(pbchIndices, rxsf, hestsf);
    [~, ~, nfmod4, mib, CellRefP] = ltePBCHDecode(rmc, pbchRx, pbchHest, nest);

    if CellRefP == 0
        fprintf("No PBCH detected for frame %d\n", frame + 1);
        continue;
    end

    rmc.CellRefP = CellRefP;
    rmc = lteMIB(mib, rmc);
    rmc.NFrame = rmc.NFrame + nfmod4;
    recFrames(frame+1) = rmc.NFrame;
    fprintf("Frame %d: PBCH decoded successfully, Frame #: %d\n", frame + 1, rmc.NFrame);
    
    rxdata = []; 

        for sf = 0:9
        if sf == 5
            continue;
        end

        rmc.NSubframe = sf;
        rxsf = rxGrid(:, frame * LFrame + sf * Lsf + (1:Lsf), :);

        % Channel estimation for current subframe
        [hestsf, nestsf] = lteDLChannelEstimate(rmc, cec, rxsf);

        % Decode PCFICH
        pcfichIndices = ltePCFICHIndices(rmc);
        [pcfichRx, pcfichHest] = lteExtractResources(pcfichIndices, rxsf, hestsf);
        [cfiBits, ~] = ltePCFICHDecode(rmc, pcfichRx, pcfichHest, nestsf);
        rmc.CFI = lteCFIDecode(cfiBits);

        % Get PDSCH symbols
        [pdschIndices, pdschIndicesInfo] = ltePDSCHIndices(rmc, rmc.PDSCH, rmc.PDSCH.PRBSet);
        [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsf, hestsf);

        % Decode PDSCH
        [rxEncodedBits, rxEncodedSymb] = ltePDSCHDecode(rmc, rmc.PDSCH, pdschRx, pdschHest, nestsf);

        % DL-SCH decode
        outLen = rmc.PDSCH.TrBlkSizes(sf + 1);
        [decbits{sf+1}, blkcrc(sf+1)] = lteDLSCHDecode(rmc, rmc.PDSCH, outLen, rxEncodedBits);

        if blkcrc(sf+1) == 0
            rxdata = [rxdata; decbits{sf+1}{:}];  % Only add if CRC passed
        else
            fprintf("CRC failed for frame %d subframe %d – skipped\n", frame+1, sf);
        end

        % Optional: collect symbols for constellation
        if sf == 0
            rxSymbols = rxEncodedSymb{1};
            % Re-encode for comparison
            txRecode = lteDLSCH(rmc, rmc.PDSCH, pdschIndicesInfo.G, decbits{sf + 1});
            txRemod = ltePDSCH(rmc, rmc.PDSCH, txRecode);
            [~, refSymbols] = ltePDSCHDecode(rmc, rmc.PDSCH, txRemod);
            txSymbols = refSymbols{1};
        end
    end

        % Concatenate decoded bits (excluding subframe 5)
    rxdata = [];
    for i = 1:length(decbits)
        if i ~= 6
            rxdata = [rxdata; decbits{i}{:}];
        end
    end

    rxDataFrame(:, frame + 1) = rxdata;
end

% Reconstruct full received bitstream
decodedRxDataStream = reshape(rxDataFrame, [], 1);

% Optional BER calculation
bitLen = min(length(decodedRxDataStream), length(trData));
bitErrorRate = comm.ErrorRate;
err = bitErrorRate(decodedRxDataStream(1:bitLen), trData(1:bitLen));

fprintf("\nBER: %.5f | Bit Errors: %d | Total Bits Compared: %d\n", err(1), err(2), err(3));

% Round bit length to nearest multiple of 8
bitLen = bitLen - mod(bitLen, 8);
bitChars = reshape(sprintf('%d', decodedRxDataStream(1:bitLen)), 8, []).';

% Convert 8-bit binary to characters
decodedText = char(bin2dec(bitChars))';

% Show result
fprintf("\nDecoded Text:\n%s\n", decodedText);

% --- Final metrics ---
signalPower = bandpower(rxWaveform);
actualSNRdB = 10 * log10(signalPower / noisePower);

bitLen = min(length(decodedRxDataStream), length(trData));
bitErrorRate = comm.ErrorRate;
err = bitErrorRate(decodedRxDataStream(1:bitLen), trData(1:bitLen));


evmRMS = sqrt(mean(abs(txSymbols - rxSymbols).^2));
evmPercent = evmRMS * 100;

drift = mean(abs(rxSymbols));

fprintf("\n--- Final Report ---\n");
fprintf("AWGN Power: %.6f\n", noisePower);
fprintf("Measured SNR: %.2f dB\n", actualSNRdB);
fprintf("BER: %.5f | Bit Errors: %d | Total Bits Compared: %d\n", err(1), err(2), err(3));
fprintf("EVM: %.2f %% RMS\n", evmPercent);
fprintf("Constellation center drift (avg): %.3f\n", drift);


% Chose output filePath
% [outputName, outputPath] = uiputfile('output.txt', 'Save Output Text As');
% if isequal(outputName,0)
%     error('Output file save canceled.');
% end
% outputFile = fullfile(outputPath, outputName);
% 
% fid = fopen(outputFile, 'w');
% fwrite(fid, decodedText);
% fclose(fid);
% 
% fprintf("\nDecoded output saved to: %s\n", outputFile);
%

outputFile = 'C:\Users\mfb36\Desktop\output.txt';
fid = fopen(outputFile, 'w');
fwrite(fid, decodedText);
fclose(fid);

fprintf("\nDecoded output saved to: %s\n", outputFile);

% Constellation Diagram
if exist('txSymbols','var') && ~isempty(txSymbols) && ...
   exist('rxSymbols','var') && ~isempty(rxSymbols)
    figure;
    scatter(real(txSymbols), imag(txSymbols), 'bo'); hold on;
    scatter(real(rxSymbols), imag(rxSymbols), 'rx');
    legend('Transmitted Symbols', 'Received Symbols');
    title('Constellation Diagram');
    xlabel('In-Phase'); ylabel('Quadrature');
    axis equal; grid on;
else
    warning('Constellation data missing. txSymbols or rxSymbols is empty.');
end

% Transmitted vs Received Bitstream Plot
if exist('trData','var') && exist('decodedRxDataStream','var')
    comparisonLength = min(length(trData), length(decodedRxDataStream));
    txBits = trData(1:comparisonLength);
    rxBits = decodedRxDataStream(1:comparisonLength);

    figure;
    subplot(2,1,1);
    stairs(txBits, 'b', 'LineWidth', 1.2);
    title('Transmitted Bitstream');
    ylim([-0.2, 1.2]); ylabel('Bit Value');
    grid on;

    subplot(2,1,2);
    stairs(rxBits, 'r', 'LineWidth', 1.2);
    title('Received Bitstream');
    ylim([-0.2, 1.2]); xlabel('Bit Index'); ylabel('Bit Value');
    grid on;
else
    warning('Cannot plot bitstream comparison. Check trData and decodedRxDataStream.');
end

figure;
pwelch(rxWaveform, [], [], [], rmc.SamplingRate, 'centered');
title('Received Signal Power Spectral Density');
