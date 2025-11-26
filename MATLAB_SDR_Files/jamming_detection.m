% Jamming Detection and Clean Frequency Selection
% This block scans nearby frequencies to detect possible jamming (interference)
% and selects the cleanest frequency (lowest power) for transmission.

% Inputs:
% - sdrReceiver: SDR receiver object (configured with initial parameters)
% - centerFrequency: Initial center frequency (Hz)
% - captureTime: Duration of each frequency scan (seconds)

% Outputs:
% - centerFrequency: Updated clean frequency (Hz) with lowest observed power
% - powerLevels: Power values (dB) at each scanned frequency for debugging/analysis

% Define a set of nearby frequencies to scan
scanFrequencies = [ ...
    centerFrequency - 6e6, ...  % -6 MHz offset
    centerFrequency - 3e6, ...  % -3 MHz offset
    centerFrequency,      ...  % original frequency
    centerFrequency + 3e6, ...  % +3 MHz offset
    centerFrequency + 6e6  ...  % +6 MHz offset
];

% Initialize vector to hold measured power levels
powerLevels = zeros(length(scanFrequencies), 1);

% Loop through each frequency, set receiver center frequency, and measure signal power
fprintf("\n--- Jamming Detection Phase ---\n");
for k = 1:length(scanFrequencies)
    sdrReceiver.CenterFrequency = scanFrequencies(k);
    rxSnapshot = capture(sdrReceiver, captureTime, "Seconds");
    powerLevels(k) = pow2db(bandpower(rxSnapshot));  
    fprintf("Scanned Frequency: %.2f MHz | Measured Power: %.2f dB\n", ...
        scanFrequencies(k)/1e6, powerLevels(k));
end

% Select the frequency with the lowest power (assumed to be least jammed)
[~, cleanIndex] = min(powerLevels);
centerFrequency = scanFrequencies(cleanIndex);

% Update both transmitter and receiver to use the selected clean frequency
sdrTransmitter.CenterFrequency = centerFrequency;
sdrReceiver.CenterFrequency = centerFrequency;

% Output the final selected frequency
fprintf("Clean frequency selected: %.2f MHz (lowest power observed)\n", centerFrequency / 1e6);
