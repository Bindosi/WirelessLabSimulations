% Step 1: Encode text message to binary
clear all
text = 'hello world';
binaryData = dec2bin(double(text), 8); % Convert characters to ASCII and then to binary
binaryData = reshape(binaryData.', 1, []); % Convert to row vector

% Step 2: Map binary data onto symbols
symbols = 2 * (binaryData.' - '0') - 1; % BPSK modulation

% Step 3: OFDM modulation
subcarriers_total = 64;
subcarriers_used = 48;
cp_size = 4;
num_symbols = ceil(length(symbols) / subcarriers_used); % Number of OFDM symbols needed
ofdm_symbols = zeros(num_symbols, subcarriers_total); % Initialize OFDM symbols matrix

% Map symbols onto subcarriers
for i = 1:num_symbols
    start_idx = (i - 1) * subcarriers_used + 1;
    end_idx = min(i * subcarriers_used, length(symbols));
    ofdm_symbols(i, 1:subcarriers_used) = symbols(start_idx:end_idx);
end

% Step 4: Add cyclic prefix
ofdm_symbols_with_cp = [ofdm_symbols(:, end-cp_size+1:end), ofdm_symbols];

% Step 5: Simulate channel effects
channel_taps = [1, 0.5, 0.3, 0.2]; % Example channel impulse response
received_signal = conv2(ofdm_symbols_with_cp, channel_taps, 'valid');

% Step 6: Add noise
SNR_dB = 10; % Signal-to-Noise Ratio in dB
noise_power = 10^(-SNR_dB / 10);
received_signal = received_signal + sqrt(noise_power) * randn(size(received_signal));

% Step 7: OFDM demodulation
received_signal_length = size(received_signal, 1);
received_ofdm_symbols = received_signal(:, cp_size+1:end);

% Step 8: Decode symbols back to text
received_symbols = reshape(received_ofdm_symbols.', 1, []);
received_binary_data = (received_symbols > 0);
received_binary_data = received_binary_data(1:length(binaryData)); % Trim extra bits
received_text = char(bin2dec(reshape(char(received_binary_data + '0'), 8, []).'));

disp(['Transmitted Text: ', text]);
disp(['Received Text: ', received_text]);
