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
W                    = 2;

% Message Source and Bits
Message = ['We hold these truths 2 be self-evident,' ...
    ' that all men are created equal, that they are ' ...
    'endowed by their Creator with certain unalienable ' ...
    'Rights, that among these are Life, Liberty and ' ...
    'the pursuit of Happiness.'];

MessageBits = reshape((dec2bin(Message,8) - '0').',1,[]);

%MessageBits = char2bin(Message);

%UsedSubcarrierIndx      = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
UsedSubcarrierIndx      = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
availableCarriers           = UsedSubcarrierIndx;
pilotCarriers               = availableCarriers(2:3:end);
availableCarriers(2:3:end)  = [];
effectiveDataCarriers   = availableCarriers;


reshapedModulatedSymbols    = reshape(MessageBits,length(effectiveDataCarriers),length(MessageBits)/length(effectiveDataCarriers));
NumSymbols                  = length(reshapedModulatedSymbols(:))/length(effectiveDataCarriers);
% Creating IFFT matrix for OFDM Modulation
to_ifft = zeros(FFTsize,NumSymbols);
% Modulating message bits to BPSK
modulatedSymbols = qammod(reshapedModulatedSymbols,2);

% Generating Pilot Bits
pilotBits              = repmat([1 0],1,length(pilotCarriers)/2);
% Modulating pilot bits and reshapinig them into pilotcarriers x Number of symbols
pilotBitsModulated     = pilotBits*2-1;
pilotSymbols           = repmat(pilotBitsModulated(:),1,NumSymbols);

% inserting pilot bits
to_ifft(pilotCarriers,:)            = pilotSymbols;
% inserting message bitss
to_ifft(effectiveDataCarriers,:)    = modulatedSymbols;

% Performing OFDM modulation IFFT
    timeDomainOFDMSymbols           = ifft(to_ifft).*sqrt(FFTsize);
% Creating our CP   
    CP = timeDomainOFDMSymbols(end-CPsize+1:end,:);
% Adding CP to the OFDM symbol
    timeDOFDMWithCp                 = [CP; timeDomainOFDMSymbols];
    timeDomainOFDMSymbolsSerial     = timeDOFDMWithCp(:).';

figure()
[Pxx_tx,F_tx] = pwelch(timeDomainOFDMSymbolsSerial,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))

% Frame Generation and filtering
Pr          = (mseq(2,6))';
Preambles   = [Pr Pr];

Frame       = [zeros(1,50) Preambles timeDomainOFDMSymbolsSerial zeros(1,50) ];
temp        = repmat(Frame,OSR,1);
TxFrame     = reshape(temp,1,length(Frame)*OSR);

% Channel
[RxSignal,H_cir] = Channeling(TxFrame.',N_taps,fs);

figure()
waterfall(abs(H_cir))
title('Plotting the Channel');
plot(abs(RxSignal))

% Synchronization
 [SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);
 
 figure()
 scatterplot(SyncPreamble)

 %Preamble
 estimated_channel      = mean(SyncPreamble.'./Preambles);
 equalizedPreambles     = SyncPreamble*conj(estimated_channel);
 scatterplot(equalizedPreambles)

 %Data
 syncOFDMmatrix         = reshape(SyncOFDMsignal,FFTsize+CPsize,[]);
 syncOFDMCPremoved      = syncOFDMmatrix(CPsize+1:length(syncOFDMmatrix),:);
 
 fdomainSymbol          = fft(syncOFDMCPremoved)/sqrt(FFTsize);
 receivedPilots         = fdomainSymbol(pilotCarriers,:);
 estimatedOFDMChannel   = receivedPilots./pilotSymbols;

% Channel Interpolation using interp
interPolated            = interp1(pilotCarriers,estimatedOFDMChannel,1:1:length(fdomainSymbol),"spline");
interPolatedChannel     = interPolated(availableCarriers,:);

waterfall(abs(interPolatedChannel));

receivedDataOFDM        = fdomainSymbol(availableCarriers,:);
equalized_frame         = receivedDataOFDM(:)./interPolatedChannel(:); 
receivedBits            = qamdemod(equalized_frame,2);

scatterplot(equalized_frame)

%Binary to character conversion
receivedMessage         = bin2char(receivedBits);

msgbox(receivedMessage,"The received Message")

mean(receivedBits(:)==MessageBits.')



