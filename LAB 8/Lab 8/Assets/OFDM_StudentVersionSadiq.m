clear
%%
% Channel Parameters
Max_Excess_Delay   = 3.2*1e-5; % in seconds

% OFDM Signal Specs
FFTsize              = 64;
df                   = 15e3; % subcarrier spacing
fs                   = df*FFTsize;
Ts                   = 1/fs;
OSR                  = 8;
NumActiveSubcarriers = 48;
% NumActiveSubcarriers32 = 32;
% NumActiveSubcarriers16 = 16;
NumSymbols           = 25;
N_taps               = floor(Max_Excess_Delay/Ts)+1;
CPsizeOpt            = ceil(Max_Excess_Delay/Ts/OSR); % Minimum CP size required (effect of OSR taken into acount)
CPsize               = CPsizeOpt;

CPsize2 = 16; % the ratio of the additional guard period

frame_RC=[]; x0 = zeros(1, NumSymbols +CPsize);
%------Raised-cosine (comparison)-------------
RRC = 0.5 + 0.5.*cos(pi + pi.*[0:CPsize2-1]./CPsize2);
LRC = 0.5 + 0.5.*cos(pi.*[0:CPsize2-1]./CPsize2);


% Message Source and Bits
Message = ['We hold these truths 2 be self-evident,' ...
    ' that all men are created equal, that they are ' ...
    'endowed by their Creator with certain unalienable ' ...
    'Rights, that among these are Life, Liberty and ' ...
    'the pursuit of Happiness'];
MessageBits = reshape((dec2bin(Message,8)-'0').',1,[]);
% MessageBinChar = dec2bin(double(char(Message)));
% [Sz1,Sz2]      = size(MessageBinChar);
% MessageBinStrn = reshape(MessageBinChar,1,Sz1*Sz2);
% MessageBits    = zeros(1,Sz1*Sz2);
% for b=1:Sz1*Sz2
%     MessageBits(b) = str2double(MessageBinStrn(b));
% end
%% Pilot bit generation and insertiom
ln = length(MessageBits);
MessagBitsWithPilot = zeros(1, ln + floor(ln / 3));

MessagBitsWithPilot(1) = MessageBits(1);
cb=1;
j = 2;
for i = 2:3:length(MessagBitsWithPilot)
    MessagBitsWithPilot(i) = cb;
    MessagBitsWithPilot(i+1:i+2) = MessageBits(j:j+1);

    cb = 1 - cb;
    j = j + 2;
end
pilotBitsCheck = MessagBitsWithPilot(2:3:length(MessagBitsWithPilot))

%%
% BPSK mappıng
BPSK_symbols = qammod(MessagBitsWithPilot,2);

BPSK_symbols_pad = [BPSK_symbols zeros(1,NumActiveSubcarriers - mod(length(BPSK_symbols),NumActiveSubcarriers))];
% BPSK_symbols_pad32 = [BPSK_symbols zeros(1,NumActiveSubcarriers32 - mod(length(BPSK_symbols),NumActiveSubcarriers32))];
% BPSK_symbols_pad16 = [BPSK_symbols zeros(1,NumActiveSubcarriers16 - mod(length(BPSK_symbols),NumActiveSubcarriers16))];
% OFDM modulation
clippingThreshold = 0.9; % Adjust as needed

UsedSubcarrierIndx = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
% UsedSubcarrierIndx32 = [FFTsize-(NumActiveSubcarriers32/2:-1:1)+1, 1:NumActiveSubcarriers32/2];
% UsedSubcarrierIndx16 = [FFTsize-(NumActiveSubcarriers16/2:-1:1)+1, 1:NumActiveSubcarriers16/2];
frame = [];
frame_clipped  = [];
% frame32 = [];
% frame16 = [];
for ii = 1:NumActiveSubcarriers:length(BPSK_symbols_pad)
    current_BPSK_symbols = BPSK_symbols_pad(ii:ii+NumActiveSubcarriers-1);
    % Clip the BPSK symbols
    current_BPSK_symbols_clipped = min(current_BPSK_symbols, clippingThreshold) + max(current_BPSK_symbols, -clippingThreshold);

    to_ifft = zeros(1,FFTsize);
    to_ifft(UsedSubcarrierIndx) = current_BPSK_symbols;
    to_ifft_clipped(UsedSubcarrierIndx) = current_BPSK_symbols_clipped;

    ofdm_symbol = ifft(to_ifft);
    ofdm_symbol_clipped = ifft(to_ifft_clipped);

    ofdm_symbol_with_cp = [ofdm_symbol(end-CPsize+1:end) ofdm_symbol];
    ofdm_symbol_with_cp_clipped = [ofdm_symbol_clipped(end-CPsize+1:end) ofdm_symbol_clipped];
    frame = [frame ofdm_symbol_with_cp];
    frame_clipped = [frame_clipped ofdm_symbol_with_cp_clipped];
    %------Raised-cosine (comparison)-------------
    Ccp = LRC.*x0(CPsize+[1:CPsize2]) + RRC.*frame(end-CPsize-CPsize2+[1:CPsize2]);
    frame_RC = [frame_RC, Ccp, frame];
    x0 = frame;


end


q=frame';


PAPR = max(abs(frame).^2) / mean(abs(frame).^2); % for 48 subC is 22.2056, for 32 subC is 15.3772, and for 16 subC is 8.9612
figure;
[Pxx, f] = pwelch(frame,[],[],[],fs,'centered', "power");
plot(f, 20*log10(Pxx));
title('Transmitted PSD');

qq=frame';
ccdf = comm.CCDF(...
    'AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true);
[ccdfy,ccdfx,avg,peak] = ccdf(qq);
figure
plot(ccdf)
legend('48 active subcarriers')

figure

[Pxx1,f]=pwelch(frame,[],[],2048,1,'twosided','centered');
plot(f,10*log10(Pxx1))
legend('Without windowing')
figure
[Pxx2,f]=pwelch(frame_clipped,[],[],2048,1,'twosided','centered');
plot(f,10*log10(Pxx2),'r--')
legend('With windowing')

ccdf1 = comm.CCDF(...
    'AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true);
[ccdfy1,ccdfx1,avg,peak] = ccdf1(frame_clipped');
figure
plot(ccdf)
hold on
plot(ccdf1)
legend('Before clipping', 'After clipping')
%% Preamble
v = mseq (2 , 7) ;
p = [ v ; v ]; % to change the matrix dimension
p = p / rms ( p ) ; % normalized the preamble

guard = zeros (50 ,1) ;
% seq =[ p ; symbs_tx ; guard ];
signal =[ guard ; p ; frame' ; guard ];


flt_stp.sps=8;
flt_stp.span=12;
flt_stp.alpha=0.1;
flt_stp . BT = 0.5;
fltr = rectpulse(signal,OSR);
fltr_up = upsample( fltr , flt_stp.sps ) ;
fltr_signal = conv ( fltr , fltr_up ) ;

figure
plot(20*log10(abs(fltr_signal)))
xlabel('Time (s)');
ylabel('Power (dB)');
title('Flitered Tx Frame')

%% Channel
% [RxSignal,H_cir] = Channeling(TxFrame,N_taps,fs);
[rx_signal,CIR]=Channeling(fltr,N_taps,fs);
RxSignal=rx_signal;
figure
plot(20*log10(abs(rx_signal)))
xlabel('Time (s)');
ylabel('Power (dB)');
title('Flitered Tx Frame')

% Plot CIR using waterfall
figure;
waterfall(abs(CIR));
xlabel('Tap Index');
ylabel('Time Index');
zlabel('Magnitude');
title('Channel Impulse Response (CIR)');

% Plot RX Frame
figure;
[CIR,rx_signal] = meshgrid(-3:.125:3);
Z = peaks(CIR,rx_signal);
waterfall(CIR,rx_signal,Z)
xlabel('Sample Index');
ylabel('Time Index');
zlabel('Magnitude');
title('Received (RX) Frame');


%% Synchronization
[SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);
%% Channel Estimation and Equalization for Preamble
% Assuming 'SyncPreamble' contains the synchronized preamble symbols

% Extract preamble symbols
preamble_symbols = SyncPreamble(1:size(126, 1), : );

% Estimate the average channel response from preamble symbols
average_channel_response = mean(preamble_symbols, 2);

% Equalize the synchronized preamble symbols using the estimated channel response
equalized_preamble_symbols = preamble_symbols ./ average_channel_response;

% Plot constellation of the equalized preamble symbols
scatterplot(equalized_preamble_symbols);
scatterplot(real(v))

fftsignal=fft(SyncOFDMsignal,48,64);

%% Channel Estimation Using Pilots
% Assuming 'OFDM_Pilot_channel' contains the received pilot symbols after synchronization

% Extract the pilot symbols from the received signal 
pilot_symbols = SyncOFDMsignal(2:3:length(SyncOFDMsignal));
pilotsOriginal = repmat([1 0],1,568/2);
pilotsOriginal(end) = [];
modulatedPilots = qammod(pilotsOriginal,2)
% Assuming 'Pilot_vector' contains the original transmitted pilot symbols

% Perform channel estimation by dividing received pilot symbols by transmitted pilot symbols
channel_estimation = pilot_symbols.' ./ modulatedPilots;

%% Interpolation for Data Carriers
pilot_index = 2:3:length(SyncOFDMsignal);
%% Interpolation for Data Carriers
% Check if there are enough pilot symbols for interpolation
if length(pilot_index) < 2
    error('Insufficient pilot symbols for interpolation. At least two distinct pilot symbols are required.');
end

% Interpolate the channel response for data carriers using interp1 function
data_carrier_indices = setdiff(1:length(SyncOFDMsignal), pilot_index); % Indices of data carriers
data_carrier_channel_response = interp1(pilot_index, channel_estimation, 1:1:length(SyncOFDMsignal), 'spline', 'extrap');

datacarrierChannel = data_carrier_channel_response(data_carrier_indices)
% Plot the channel frequency response for data carriers
figure;
plot(abs(data_carrier_channel_response), 'b', 'LineWidth', 2);
title('Channel Frequency Response for Data Carriers');
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;


figure
plot(20*log10(abs(data_carrier_channel_response)))
xlabel('Time (s)');
ylabel('Power (dB)');
title('Channel Frequency Response for Data Carriers')
% OFDM demodulation, Channel estimation, equalization and detection


equalized_symbols =SyncOFDMsignal(data_carrier_indices)./ datacarrierChannel.';


scatterplot(equalized_symbols(:)); % Convert equalized_symbols to a vector
title('Constellation of Equalized Symbols');


scatterplot(equalized_preamble_symbols(:)); % Convert equalized_preamble_symbols to a vector
title('Constellation of Equalized Preamble Symbols');

% Assuming received_symbols_freq contains the demodulated received symbols

% Demodulate the received symbols to obtain the bits
received_bits = real(equalized_symbols) > 0; % Assuming BPSK modulation

receivedMessage = char(reshape(bin2dec(reshape(char(received_bits(1:1128) +'0'),8,[]).'),1,[]))


% % Convert the bits into characters
% num_chars = floor(length(received_bits) / 8); % Number of characters
% received_bits_truncated = received_bits(1:num_chars * 8); % Truncate to multiples of 8
% received_bits_reshaped = reshape(received_bits_truncated, 8, []).'; % Reshape into bytes
% received_bytes = bin2dec(received_bits_reshaped); % Convert bytes from binary to decimal
% received_message = char(received_bytes.'); % Convert decimal values to characters

% Plot Bit Error Rate (BER)
SNRdB = 0:2:20; % Signal-to-Noise Ratio (SNR) in dB
BER = zeros(size(SNRdB));

for i = 1:length(SNRdB)
    % Add noise to the received bits based on the SNR
    received_bits_noisy = awgn(SyncOFDMsignal, SNRdB(i), 'measured');
    
    % Calculate the number of errors
    num_errors = sum(BPSK_symbols ~= received_bits_noisy);
    
    % Calculate the Bit Error Rate (BER)
    BER(i) = num_errors/ length(BPSK_symbols);
end

% Plot BER vs. SNR
figure;
semilogy(SNRdB, BER, '-o');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate (BER) vs. SNR');