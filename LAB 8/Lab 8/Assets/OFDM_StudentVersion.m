clear; clc;

% Channel Parameters
Max_Excess_Delay   = 3.2*1e-5; % in seconds

% OFDM Signal Specs
FFTsize              = 64;
df                   = 15e3; % subcarrier spacing
fs                   = df*FFTsize;
Ts                   = 1/fs;
OSR                  = 8;
NumActiveSubcarriers = 48;
NumSymbols           = 25;
N_taps               = floor(Max_Excess_Delay/Ts)+1;
CPsizeOpt            = ceil(Max_Excess_Delay/Ts/OSR); % Minimum CP size required (effect of OSR taken into acount)
CPsize               = CPsizeOpt;

% Message Source and Bits
Message = ['We hold these truths 2 be self-evident,' ...
    ' that all men are created equal, that they are ' ...
    'endowed by their Creator with certain unalienable ' ...
    'Rights, that among these are Life, Liberty and ' ...
    'the pursuit of Happiness'];

MessageBinChar = dec2bin(double(char(Message)));
[Sz1,Sz2]      = size(MessageBinChar);
MessageBinStrn = reshape(MessageBinChar,1,Sz1*Sz2);
MessageBits    = zeros(1,Sz1*Sz2);
for b=1:Sz1*Sz2
    MessageBits(b) = str2double(MessageBinStrn(b));
end

pilotBitsIndices = [2];
for k=5:3:length(MessageBits)
    pilotBitsIndices = [pilotBitsIndices;k];
end

totalLength = length(MessageBits)+ length(pilotBitsIndices);

% Pilot bit generation and insertiom
PilotBits = [1;0];

for k = 1:length(pilotBitsIndices)/2
    PilotBits = [PilotBits; 1;0];
end

PilotBits(end) = [];
MessageAndPilots = zeros(totalLength,1);

% inserting the pilots into the message bits
MessageAndPilots(pilotBitsIndices) = PilotBits;

% Here we are getting the index positions for the message bits
preMessageIndices                  = 1:totalLength; % generating all indexes
preMessageIndices(pilotBitsIndices) = [];       % removing all indexes used for pilots
MessageAndPilots(preMessageIndices) = MessageBits;  % Inserting the data bits into the index locations

% Modulating the generated bits into BPSK symbols
BPSK_modulated_symbols =  pskmod(MessageAndPilots,2)';

% OFDM modulation
 UsedSubcarrierIndx     = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
 AllSubCarriers         = 1:FFTsize;
 AllSubCarriers(UsedSubcarrierIndx) = [];
 emptyCarriers          =  AllSubCarriers;
 ofdmSymbolFFTsize      =  zeros(FFTsize,1); 

 
 % X_blocks = reshape(X_padded,nfft,length(X_padded)/nfft);
 num_symbols = ceil(length(BPSK_modulated_symbols) / NumActiveSubcarriers); % Number of OFDM symbols needed
 ofdm_symbols = zeros(num_symbols, FFTsize); % Initialize OFDM symbols matrix
 
 % Map symbols onto subcarriers
for i = 1:num_symbols
    start_idx = (i - 1) * UsedSubcarrierIndx + 1;
    end_idx = min(i * UsedSubcarrierIndx, length(BPSK_modulated_symbols));
    ofdm_symbols(i, 1:UsedSubcarrierIndx) = BPSK_modulated_symbols(start_idx:end_idx);
end

% Step 4: Add cyclic prefix
ofdm_symbols_with_cp = [ofdm_symbols(:, end-CPsize+1:end), ofdm_symbols];


% Frame Generation and filtering
ifft_OFDM_Symbol_CP_Oversampled = ofdm_symbols_with_cp * OSR;
TxFrame =       ifft_OFDM_Symbol_CP_Oversampled;


% % Channel
 [RxSignal,H_cir] = Channeling(TxFrame,N_taps,fs);


% % Synchronization
 [SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);


% OFDM demodulation, Channel estimation, equalization and detection

