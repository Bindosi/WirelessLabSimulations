
mod_stp.T = 'QAM';
mod_stp.M = 2;
mod_stp.N = 512;


flt_stp.T = 'RRC';
flt_stp.sps = 8;
flt_stp.span = 12;
flt_stp.alpha = 0.9;
flt_stp.BT = 0.7;
flt_stp.mtchd = true;

%Hardware Setup
hw_stp.Tx_G  = -10;
hw_stp.Rx_G  = 20;
hw_stp.F     = 2.35e9;
hw_stp.R     = 1e6;
hw_stp.Tx_ID = 'sn:104473541196000cf7ff14008dd81e43f6';
%hw_stp.Tx_ID = 'sn:10447354119600060d001800cf281e583b';
hw_stp.Rx_ID ='sn:104473dc599300131200210082672a4170';
hw_stp.N = 30e3;

% setting up transmit object
        Tx  = sdrtx('Pluto', 'Gain',hw_stp.Tx_G,...
            'CenterFrequency',hw_stp.F,...
            'BasebandSampleRate',hw_stp.R,...
            'RadioID',hw_stp.Tx_ID);

% setting up receive object
        Rx = sdrrx('Pluto', 'SamplesPerFrame', hw_stp.N,...
            'OutputDataType','double',...
            'CenterFrequency',hw_stp.F,...
            'BasebandSampleRate',hw_stp.R,...
            'GainSource','Manual',... %AGC Fast Attack
            'Gain', hw_stp.Rx_G,...
            'RadioID',hw_stp.Rx_ID);

%% question 2, transmitting a tone

t                    = 0:1/hw_stp.R:1-1/hw_stp.R;          % Time vector for 1 second
baseband_frequency  = 1000; 
baseband_frequency1  = 6000; 
% Frequency of the tone (1 kHz)
sinusoidal_tone     = sin(2*pi*baseband_frequency*t) +sin(2*pi*baseband_frequency1*t);     % Generate a sine wave tone

% Viewing the generated signal wave
figure
plot(t(1:10000), sinusoidal_tone(1:10000))

% checking our singletone signal
figure()
[Pxx_tx_q2,F_tx_q2] = pwelch(sinusoidal_tone,[],[],[],hw_stp.R,'centered');
plot(F_tx_q2,10*log10(Pxx_tx_q2));
grid on

%Transmiting a tone from SDR
Tx.transmitRepeat(sinusoidal_tone'); 
release(Tx);
%% Question 3: Measuring Received Signal Power with los at two different distances

% Receiving
[smpls_rx_q2, prmb_rx_q2, smpls_dt_q2, symbs_rx_q2] = receive(Rx,mod_stp, flt_stp,hw_stp);

release(Rx);
figure()
[Pxx_tx_q2,F_tx_q2] = pwelch(smpls_rx_q2,[],[],[],hw_stp.R,'centered');
plot(F_tx_q2,10*log10(Pxx_tx_q2));
grid on

%% question 4: Change distances with non los links between Tx and RX
figure()
[Pxx_tx_q4,F_tx_q4] = pwelch(smpls_rx_q4,[],[],[],hw_stp.R,'centered');
plot(F_tx_q4,10*log10(Pxx_tx_q4));
grid on

%% question 5: comment on the results in 3 and 4

%% question 6: Repeat 3 to 5 with another operating frequency ie: 
% compare path loss and variations at different frequencies
hw_stp.F        = 900e6;

%% FREQUENCY SELECTIVITY
%% Question 7: Setting up variables for frequency selectivity 

% gnerate preamble
Preamble_1_q7      = mseq(2,7);
Preamble_2_q7      = mseq(2,7);

% Generate Data Bits
data_q7            =  randi([0 mod_stp.M-1],1,mod_stp.N);

% Generate BPSK symbols
mod_symbols_q7     = pskmod(data_q7, mod_stp.M); % it has 0 complex part
mod_symbols_    = 2 * data_q7 - 1;            % no 0 complex part
guard_samples   = zeros(1,200); 
trasmit_frame_q7   = [guard_samples Preamble_1_q7' Preamble_2_q7' mod_symbols_q7 guard_samples];

% Upsampling
tx_frame_up__q7 = upsample(trasmit_frame__q7,flt_stp.sps);

% creatıng our fılter
filter = rcosdesign(flt_stp.alpha, flt_stp.span, flt_stp.sps);

% pefroming convolution to get the filtered frame
filteredFrame_q7 = conv(filter, tx_frame_up__q7);

% testing the look 
figure()
[Pxx_tx__q7,F_tx_q7] = pwelch(filteredFrame_q7,[],[],[],hw_stp.R,'centered');
plot(F_tx_q7,10*log10(Pxx_tx_q7));
grid on

%% Question 9: Testing the prepared samples using different sample rates

% Change hw.stp.R = 0.5Msps, 1Msps, 10Msps, 15Msps;
%[frm,symbs_tx] = transmit(Tx,mod_stp,flt_stp);
 Tx.transmitRepeat(filteredFrame_q14); 

% Receiving
[smpls_rx_q9, prmb_rx_q7, smpls_dt_q7, symbs_rx_q7] = receive(Rx,mod_stp, flt_stp,hw_stp);

figure()
[Pxx_tx_q9,F_tx_q9] = pwelch(smpls_rx_q7,[],[],[],hw_stp.R,'centered');
plot(F_tx_q9,10*log10(Pxx_tx_q9));
grid on

%% Question 10
% Change hw.stp.R = 10MSps, 20MSps, 30MSps);
%  Here we use VSG and VSA
%% TIME SELECTIVITY
% Question 11 Connect Antennas to the reverberation chambers st 2.4 Ghz transmit fx
% Question 12 We also use VSG and VSA turn on the fan inside to create Doppler spread
% Question 13 Repeat Question 11 and 2 at 900Mhz

%% CHANNEL ESTIMATION
% gnerate preamble
Preamble_1_q14      = mseq(2,7);
Preamble_2_q14      = mseq(2,7);
Preamble_q14        = [Preamble_1_q14' Preamble_2_q14'];
% Setting modulation order
mod_stp.M_q14       = 4;
% Generate Data Bits
data_q14 =  randi([0 mod_stp.M_q14-1],1, mod_stp.N);

% Generate BPSK symbols
mod_symbols_q14     = pskmod(data_q14, mod_stp.M_q14, pi/mod_stp.M_q14);
scatterplot(mod_symbols_q14)
guard_samples       = zeros(1,200); 
trasmit_frame_q14   = [guard_samples Preamble_1_q14.' Preamble_2_q14.' mod_symbols_q14 guard_samples];
scatterplot(trasmit_frame_q14)
% Upsampling
tx_frame_up_q14     = upsample(trasmit_frame_q14,flt_stp.sps);

% creatıng our fılter
filter              = rcosdesign(flt_stp.alpha, flt_stp.span, flt_stp.sps);

% pefroming convolution to get the filtered frame
filteredFrame_q14   = conv(filter, tx_frame_up_q14);

% % testing the look 
figure()
[Pxx_tx_q14,F_tx_q14] = pwelch(filteredFrame_q14,[],[],[],hw_stp.R,'centered');
plot(F_tx_q14,10*log10(Pxx_tx_q14));
grid on

%Transmitting the signal 
 Tx.transmitRepeat(filteredFrame_q14); 

% Receiving
[smpls_rx_q14, prmb_rx_q14, smpls_dt_q14, symbs_rx_q14] = receive(Rx,mod_stp, flt_stp,hw_stp);
% % smpls_rx  = received raw data samples,
% % prmb_rx   = received preamble samples,
% % smpls_dt  = received processed data samples,
% % symbs_rx  = received data samples
% 
release(Tx);
release(Rx);

% channel             = rand(5,1)+1j*rand(5,1);
% received_samples    = conv(channel, filteredFrame_q14);
% 
% down                =downsample(received_samples(52+length(guard_samples)*8:length(received_samples)-length(guard_samples)*8-48-1),8)
% scatterplot(down)
% prmb_rx_q14          =downsample(received_samples(52+length(guard_samples)*8:length(received_samples)-8*length(data_q14)-length(guard_samples)*8-49),8)
% scatterplot(prmb_rx_q14)
% downsamp            =downsample(received_samples(8*length(Preamble_q14)+52+length(guard_samples)*8:length(received_samples)-length(guard_samples)*8-48-1),8)
% scatterplot(downsamp)
% Channel Estimation
prmb_rx_q14 = downforpre;

estimated_h         = mean(downforpre/Preamble_q14);       
estimatedChannel    = conv(prmb_rx_q14, fliplr(Preamble_q14)); % Convolution Approach
channel_check       = isequal(estimated_h,channel);

% Equalization of received symbols to get the original symbols
eq_smpls_rx         = smpls_rx_q14 * conj(estimated_h);
scatterplot(eq_smpls_rx)
%1
firstmet            = pskdemod(eq_smpls_rx,4,pi/4);
scatterplot(firstmet)

% Creating reference constellation points for comparison
table               = pskmod([0 1 2 3], mod_stp.M_q14, pi/mod_stp.M_q14);
table               = table ./ rms(table);


% Detecting the symbols using Maximum Likelihood detection
% Define QPSK constellation points
constellation       = pskmod(0:3,4,pi/4); % QPSK constellation

% Initialize variables
ml_symbols          = zeros(1, length(eq_smpls_rx));
detected_data       = zeros(1, length(eq_smpls_rx));

% Perform ML detection
for i               = 1:length(received_symbols)
    % Calculate Euclidean distance to each constellation point
    distances       = abs(received_symbols(i) - table);
    
    % Determine the index of the nearest constellation point (ML detection)
    [~, idx]        = min(distances);
    
    % Map index back to QPSK symbol
    ml_symbols(i)   = table(idx);

    % Map detected data to original data
    detected_data(i)    = mod(idx-1, M);
end
scatterplot(detected_data)
scatterplot(ml_symbols)



