% Initialize transmitter (Tx) and receiver (Rx) objects
Tx = sdrtx('Pluto');
Rx = sdrrx('Pluto');

% PATH-LOSS EXPONENT
% Specify different distances (in meters)
distances = [10, 20, 30]; 

% Specify different frequencies (in Hz)
frequencies = [2.35e9, 2.4e9, 2.45e9];

fprintf('PATH-LOSS EXPONENT\n');
for dist = distances
    for freq = frequencies
        % Set parameters for Tx and Rx objects
        Tx.CenterFrequency = freq;
        Tx.BasebandSampleRate = 1e6;
        Tx.RadioID = findPlutoRadio;
        
        Rx.CenterFrequency = freq;
        Rx.BasebandSampleRate = 1e6;
        Rx.SamplesPerFrame = 30e3;
        Rx.GainSource = 'Manual';
        Rx.Gain = 20;
        Rx.OutputDataType = 'double';
        Rx.RadioID = findPlutoRadio;
        
        % Transmit signal
        transmitRepeat(Tx, rand(1e5,1)); % Random tone signal with 1e5 samples
        
        % Capture received signal
        rxSignal = Rx();
        
        % Measure received signal power
        receivedPower = mean(abs(rxSignal).^2); % Average power
        
        % Record and compare received signal power for this distance and frequency
        fprintf('Distance: %d meters, Frequency: %.2f GHz, Received Power: %.2f dBm\n', dist, freq/1e9, 10*log10(receivedPower));
    end
end

% FREQUENCY-SELECTIVITY
% Specify sample rates (in samples/s)
sampleRates = [0.5e6, 1e6, 10e6, 15e6];

fprintf('\nFREQUENCY-SELECTIVITY\n');
for rate = sampleRates
    % Set sample rate for Tx object
    Tx.BasebandSampleRate = rate;
    
    % Generate preamble symbols using m-sequence method
    m = mseq(2, 7, 1);
    preamble = [m; m];

    % Generate data symbols
    data = randi([0 1], 512, 1); % For BPSK, we generate binary data

    % Modulate data symbols using BPSK modulation
    symb = pskmod(data, 2); % BPSK modulation with 2 symbols

    % Create guard samples
    guard_samples = zeros(200, 1);

    % Concatenate preamble, data symbols, and guard samples
    frame = [guard_samples; preamble; symb; guard_samples];

    % Filter the frame samples using RRC filter
    roll_off_factor = 0.1;
    oversampling_ratio = 8;
    RRC_filter = rcosdesign(roll_off_factor, 6 * oversampling_ratio, oversampling_ratio, 'sqrt');
    filtered_frame = filter(RRC_filter, 1, frame);

    % Normalize filtered frame to unit energy
    filtered_frame = filtered_frame / norm(filtered_frame);

    % Transmit frame
    transmitRepeat(Tx, filtered_frame);
    
    % Capture signal
    rxSignal = receive(Rx);
    
    % Plot power spectral density
    [psd, freq] = periodogram(rxSignal, rectwin(length(rxSignal)), [], Rx.BasebandSampleRate, 'centered');
    plot(freq, 10*log10(psd));
    title('Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
    pause(1); % Pause to allow time for plot to render
end

% TIME-SELECTIVITY
% Connect antennas to the chamber and devices (Assuming hardware setup)

fprintf('\nTIME-SELECTIVITY\n');
% Turn on fan inside the chamber to create Doppler spread
% (Assuming control of fan and chamber setup)

% CHANNEL ESTIMATION
fprintf('\nCHANNEL ESTIMATION\n');
% Assuming preamble symbols and data symbols are available
preambleSymbols = randi([0,1], 1, 2e7-1); % Two m-sequences
dataSymbols = randi([0,1], 1, 512); % 512 BPSK symbols

% Transmit frame with QPSK modulation
transmitRepeat(Tx, qpskModulate([preambleSymbols, dataSymbols]));

% Capture received preamble
rxPreamble = receive(Rx);

% Estimate channel using preamble
channelEstimate = estimateChannel(rxPreamble, preambleSymbols);

% Compensate for channel effects
rxData = rxData ./ channelEstimate;

% ML symbol recovery
recoveredSymbols = mlDetector(rxDataCompensated);

% Clean up
release(Tx);
release(Rx);
