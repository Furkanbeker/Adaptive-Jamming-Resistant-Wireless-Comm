% Adaptive Modulation Selection
% This block selects the most suitable modulation scheme (QPSK, 16QAM, or 64QAM)
% based on the estimated SNR from the received LTE downlink signal.

% Inputs:
% - hest: Channel estimation matrix (from lteDLChannelEstimate)
% - noisePower: Known or measured noise power
% - rmc: LTE RMC configuration object (from lteRMCDL)

% Output:
% - Updated 'rmc' structure with the selected modulation format (QPSK / 16QAM / 64QAM).
% - Displayed log of selected modulation and corresponding estimated SNR value.

% Estimate average channel power from channel estimation
avgChannelPower = mean(abs(hest(:)).^2);  
estimatedSNRdB = 10 * log10(avgChannelPower / noisePower);  

% Select modulation based on SNR thresholds
if estimatedSNRdB < 8
    selectedModulation = 'QPSK';
elseif estimatedSNRdB < 15
    selectedModulation = '16QAM';
else
    selectedModulation = '64QAM';
end

% Update the LTE configuration with the selected modulation
rmc.PDSCH.Modulation = {selectedModulation};

% Output the selected modulation scheme and SNR
fprintf("Adaptive modulation selected: %s (Estimated SNR: %.2f dB)\n", selectedModulation, estimatedSNRdB);
