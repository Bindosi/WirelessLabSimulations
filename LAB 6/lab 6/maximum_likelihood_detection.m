

% Define channel parameters
SNR_dB = 5; % Signal-to-Noise Ratio (SNR) in dB
SNR = 10^(SNR_dB/10); % Convert SNR to linear scale



% Generate random binary data (2 bits per symbol)
data_bits = randi([0, 1], 1, 512); % 512 symbols with 2 bits each

% Map binary data to QPSK symbols
qpsk_symbols = pskmod(data_bits,4,pi/4); % QPSK modulation

% Add AWGN to QPSK symbols
received_symbols = awgn(qpsk_symbols, SNR);

scatterplot(received_symbols)
% Define QPSK constellation points
constellation = pskmod(0:3,4,pi/4); % QPSK constellation

% Initialize variables
ml_symbols = zeros(1, length(received_symbols));

% Perform ML detection
for i = 1:length(received_symbols)
    % Calculate Euclidean distance to each constellation point
    distances = abs(received_symbols(i) - constellation);
    
    % Determine the index of the nearest constellation point (ML detection)
    [~, idx] = min(distances);
    
    % Map index back to QPSK symbol
    ml_symbols(i) = constellation(idx);
end


% Calculate Bit Error Rate (BER)
num_errors = sum(ml_symbols ~= qpsk_symbols);
ber = num_errors / (512 * 2); % Total number of bits transmitted = 512 * 2

fprintf('Bit Error Rate (BER): %f\n', ber);


