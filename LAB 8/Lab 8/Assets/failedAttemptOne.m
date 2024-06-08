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

MessageBitsIndices  =   1:length(MessageBits);

% Determining Pilot bits indices
PilotBitIndices     =   MessageBitsIndices(2:3:end);

% lenght of the pilot bits
pilotLength         =   length(PilotBitIndices);

% legth of data and pilot bits
totalLength = pilotLength   + length(MessageBits);

pilotBitsIndices = [2]
for k=5:3:totalLength
    pilotBitsIndices = [pilotBitsIndices;k]
end


% Pilot bit generation and insertiom
PilotBits = [];

for k = 1:length(pilotBitsIndices)/2
    PilotBits = [PilotBits; 1;0];
end
MessageAndPilots = zeros(totalLength,1);

% inserting the pilots into the message bits
MessageAndPilots(pilotBitsIndices) = PilotBits;

% Here we are getting the index positions for the message bits
preMessageIndices                  = [1:totalLength]; % generating all indexes
preMessageIndices(pilotBitsIndices) = [];       % removing all indexes used for pilots
MessageAndPilots(preMessageIndices) = MessageBits;  % Inserting the data bits into the index locations

% OFDM modulation
 UsedSubcarrierIndx = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
 

 % indexes of ofdm carriers for pilot bits

% Frame Generation and filtering



% % Channel
% [RxSignal,H_cir] = Channeling(TxFrame,N_taps,fs);


% % Synchronization
% [SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);


% OFDM demodulation, Channel estimation, equalization and detection

