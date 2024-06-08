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
CPsize               = 1;
W                    = 1;

% Message Source and Bits
Message = ['We hold these truths 2 be self-evident,' ...
    ' that all men are created equal, that they are ' ...
    'endowed by their Creator with certain unalienable ' ...
    'Rights, that among these are Life, Liberty and ' ...
    'the pursuit of Happiness My Boss.'];

MessageBits = reshape((dec2bin(Message,8) - '0').',1,[]);

%MessageBits = char2bin(Message);

%UsedSubcarrierIndx      = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
UsedSubcarrierIndx      = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
availableCarriers           = UsedSubcarrierIndx;
pilotCarriers               = availableCarriers(2:3:end);
availableCarriers(2:3:end)  = [];
effectiveDataCarriers       = availableCarriers;

numberOfZerosToIncrease                                    = [];
if mod(length(MessageBits),length(effectiveDataCarriers)) ~= 0   
numberOfZerosToIncrease     = floor(length(MessageBits)/length(effectiveDataCarriers))*length(effectiveDataCarriers)... 
                              + length(effectiveDataCarriers) - length(MessageBits);
end

zeropaddingMessage          = [MessageBits zeros(1,numberOfZerosToIncrease)];
reshapedModulatedSymbols    = reshape(zeropaddingMessage,length(effectiveDataCarriers),length(zeropaddingMessage)/length(effectiveDataCarriers));

NumSymbols                  = length(reshapedModulatedSymbols(:))/length(effectiveDataCarriers);
% Creating IFFT matrix for OFDM Modulation

to_ifft = zeros(FFTsize,NumSymbols);
% Modulating message bits to BPSK
modulatedSymbols = qammod(reshapedModulatedSymbols,2);


% Generating Pilot Bits
pilotBits              = repmat([1 0],1, floor(length(pilotCarriers)/2));

 if length(pilotBits) < length(pilotCarriers)
     if pilotBits(end) == 1
        pilotBits       = [pilotBits 0];
     else
         pilotBits      = [pilotBits 1];
     end
 else 
     if length(pilotBits)>length(pilotCarriers)  
         pilotBits = pilotBits(1:length(pilotCarriers));
    end
 end

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

%% Calculating the PAPR of the signal
Pmax = max((abs(timeDomainOFDMSymbolsSerial)).^2);
Pavg = mean((abs(timeDomainOFDMSymbolsSerial)).^2);
PAPR = Pmax/Pavg;

%% Pwelch Plot for spectrum of the signal
[Pxx_tx,F_tx] = pwelch(timeDomainOFDMSymbolsSerial,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
xlabel('Frequency');
ylabel('Amplitude');
title('Frequency Spectrum of OFDM Symbols');


%% Plotting the CCDF function of the frame
ccdf = comm.CCDF('NumPoints', 1000000);
[ccdf_values, ccdf_bins] = step(ccdf, timeDomainOFDMSymbolsSerial');
figure()
semilogy(ccdf_bins, ccdf_values);
xlabel('Amplitude');
ylabel('Probability(%)');
title('Complementary Cumulative Distribution Function (CCDF)');

%% Performing windowing
G2=ceil(Max_Excess_Delay/Ts/OSR)

RRC = 0.5 + 0.5.*cos(pi + pi.*[0:G2-1]./G2);
LRC = 0.5 + 0.5.*cos(pi.*[0:G2-1]./G2);

Ccp1= RRC.'.*timeDOFDMWithCp(1+(end-CPsize)-G2:end-CPsize,1);
Ccp2= LRC.'.*timeDOFDMWithCp(1:G2,1:end-1)+RRC.'.*timeDOFDMWithCp(1+end-CPsize-G2:end-CPsize,2:end);
Ccp3= LRC.'.*timeDOFDMWithCp(1:G2,end);

Frame_RC1=[Ccp1; timeDOFDMWithCp(:,1)];
Frame_RC2=[Ccp2; timeDOFDMWithCp(:,2:end)];
Frame_RC2=Frame_RC2(:);
Frame_RC=[Frame_RC1; Frame_RC2; Ccp3]

% Pwelch plot with Windowing
figure()
[Pxx_tx,F_tx] = pwelch(Frame_RC(:),[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
hold on
% Pwelch Plot without windowing
[Pxx_tx,F_tx] = pwelch(timeDomainOFDMSymbolsSerial,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
legend('Windowed','Not Windowed')
xlabel('Frequency');
ylabel('Amplitude');
title('Frequency Spectrum of OFDM Symbols');
hold off
%% Clipping to improve PAPR
x_mag       =   abs(timeDomainOFDMSymbolsSerial);
x_max       =   0.5*max(x_mag);
x_mean      =   0.4*mean(x_mag);

for j=1:length(timeDomainOFDMSymbolsSerial)                                  %Clipping the signals above threshold(here 70% of original value)
    if(x_mag(j)>x_max)
        x_clip(j)   =   x_max;
    else
        x_clip(j)   =   x_mag(j);
    end
    if(x_clip(j)<x_mean)
        x_clip(j)   =   x_mean;
    else
        x_clip(j)   =   x_clip(j);
    end
end
PmaxClipped                = max((abs(x_clip)).^2);
PavgClipped                = mean((abs(x_clip)).^2);

paprClipped                = PmaxClipped/PavgClipped

%plotting the symbols after clipping
figure()
plot(abs(x_clip))
title('OFDM Time Domain Signal With Clipping')
ylabel('Amplitude')
xlabel('Samples')
hold on
figure()
plot(abs(timeDomainOFDMSymbolsSerial))
title('OFDM Time Domain Signal without Clipping')
ylabel('Amplitude')
xlabel('Samples')
hold off


%% Transmitter Part Frame Generation and filtering
Pr          = (mseq(2,6))';
Preambles   = [Pr Pr];

Frame       = [zeros(1,50) Preambles timeDomainOFDMSymbolsSerial zeros(1,50) ];
temp        = repmat(Frame,OSR,1);
TxFrame     = reshape(temp,1,length(Frame)*OSR);

% Plotting the Tx frame
figure()
plot(20*log10((abs(TxFrame))));
title('The Time Domain OFDM Signal')
ylabel('Amplitude')
xlabel('Samples')
% Channel
[RxSignal,H_cir] = Channeling(TxFrame.',N_taps,fs);

figure()
waterfall(abs(H_cir))
ylabel('Number of Taps')
xlabel('Symbols')
title('Plotting the Channel');

figure()
plot(abs(RxSignal))
ylabel('Amplitude')
xlabel('Samples')
title('Received Symbols');

% Synchronization
 [SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);
 
 scatterplot(SyncPreamble)

 %Preamble
 estimated_channel      = mean(SyncPreamble.'./Preambles);
 equalizedPreambles     = SyncPreamble*conj(estimated_channel);
 scatterplot(equalizedPreambles)

 %Data
 syncOFDMmatrix         = reshape(SyncOFDMsignal,FFTsize+CPsize,[]);
 syncOFDMCPremoved      = syncOFDMmatrix(CPsize+1:length(syncOFDMmatrix(:,1)),:);
 
 fdomainSymbol          = fft(syncOFDMCPremoved)/sqrt(FFTsize);
 receivedPilots         = fdomainSymbol(pilotCarriers,:);
 estimatedOFDMChannel   = receivedPilots./pilotSymbols;

% Channel Interpolation using interp
interPolated            = interp1(pilotCarriers,estimatedOFDMChannel,1:1:length(fdomainSymbol(:,1)),"spline");
interPolatedChannel     = interPolated(availableCarriers,:);

figure()
waterfall(abs(interPolatedChannel));

receivedDataOFDM        = fdomainSymbol(availableCarriers,:);
equalized_frame         = receivedDataOFDM(:)./interPolatedChannel(:); 
receivedBits            = qamdemod(equalized_frame,2);
receivedMessageBits     = receivedBits(1:length(MessageBits)); % removing zeros

scatterplot(equalized_frame)
xlim([-2,2])
ylim([-2,2])
%Binary to character conversion
receivedMessage         = bin2char(receivedMessageBits);

msgbox(receivedMessage,"The received Message")

mean(receivedMessageBits(:)==MessageBits.')



