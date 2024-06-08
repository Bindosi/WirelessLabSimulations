clear; clc;

% Channel Parameters
Max_Excess_Delay   = 3.2*1e-8; % in seconds

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
Message = ['Hello Yunus'];

MessageBits = char2bin(Message);

UsedSubcarrierIndx      = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];

availableCarriers           = UsedSubcarrierIndx;
pilotCarriers               = availableCarriers(2:3:end);
availableCarriers(2:3:end)  = [];
effectiveDataCarriers   = availableCarriers;

numberOfZerosToIncrease     = floor(length(MessageBits)/length(effectiveDataCarriers))*length(effectiveDataCarriers)... 
                              + length(effectiveDataCarriers) - length(MessageBits);

zeropaddingMessage = [MessageBits zeros(1,numberOfZerosToIncrease)];
reshapedModulatedSymbols = reshape(zeropaddingMessage,length(effectiveDataCarriers),length(zeropaddingMessage)/length(effectiveDataCarriers));

modulatedSymbols = qammod(reshapedModulatedSymbols,2);

% modulatedDataBits = zeropaddingMessage*2-1;
% modulatedDataBits = real(pskmod(zeropaddingMessage,2))

NumSymbols        = length(modulatedSymbols(:))/length(effectiveDataCarriers);

to_ifft = zeros(FFTsize,NumSymbols);

pilotBits              = repmat([1 0],1,length(pilotCarriers)/2);
pilotBitsModulated     = pilotBits*2-1;
pilotSymbols           = repmat(pilotBitsModulated(:),1,NumSymbols);


to_ifft(pilotCarriers,:) = pilotSymbols;


to_ifft(effectiveDataCarriers,:) = modulatedSymbols;


    timeDomainOFDMSymbols = ifft(to_ifft).*sqrt(FFTsize);

    CP = timeDomainOFDMSymbols(end-CPsize+1:end,:);

    timeDOFDMWithCp = [CP; timeDomainOFDMSymbols];

    timeDomainOFDMSymbolsSerial = timeDOFDMWithCp(:).';

figure()
plot(20*log(abs(timeDomainOFDMSymbolsSerial).^2))

Pmax = max((abs(timeDomainOFDMSymbolsSerial)).^2);
Pavg = mean((abs(timeDomainOFDMSymbolsSerial)).^2);

PAPR = Pmax/Pavg;
figure()
[Pxx_tx,F_tx] = pwelch(timeDomainOFDMSymbolsSerial,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))


ccdf = comm.CCDF('NumPoints', 1000000);
[ccdf_values, ccdf_bins] = step(ccdf, timeDomainOFDMSymbolsSerial');
figure()
semilogy(ccdf_bins, ccdf_values);
xlabel('Amplitude');
ylabel('Probability(%)');
title('Complementary Cumulative Distribution Function (CCDF)');
% %% Widnowing
% %approach 1 Alphan's
% LLC = 1/2*(cos(linspace(-pi,0,W+2))+1);
% plot(LLC)
% 
% LLC(1) = [];
% LLC(end) = [];
% 
% W1 = eye(FFTsize+CPsize);
% W1(1:W,1:W) = diag(LLC);
% 
% RRC = 1/2*(cos(linspace(0,pi,W+2))+1);
% plot(RRC)
% 
% RRC(1) = [];
% RRC(end) = [];
% W2 = zeros(FFTsize+CPsize);
% 
% W2(1:W,(CPsize+(1:W))) = diag(RRC);
% 
% xWOFDM = W1*timeDOFDMWithCp  + W2* timeDOFDMWithCp;
% 
% if (1) % OOB
%     xA = timeDOFDMWithCp(:);
%     Nwin = 1024;
%     Noverlap= floor(Nwin*3.2/4);
%     [pX,f] = pwelch(xA(:), Nwin, Noverlap,[],FFTsize,'twosided');
%     pX = pX * FFTsize;
%      h100 = figure(100);
%     plot(f-FFTsize/2,fftshift(10*log10(pX)) , 'r', 'displayname', 'OFDM' )
%     hold on
% 
%     xA = xWOFDM(:).';
%     Nwin = 1024;
%     Noverlap= floor(Nwin*3.2/4);
%     [pX,f] = pwelch(xA(:), Nwin, Noverlap,[],FFTsize,'twosided');
%     pX = pX * FFTsize;
%     h100 = figure(100);
%     plot(f-FFTsize/2,fftshift(10*log10(pX)) , 'g', 'displayname', 'Windowing & OFDM' )
%     hold on
% 
% 
%     grid on
%     xlabel('Subcarrier Index')
%     ylabel('Power spectral density [dBm/Hz]')
%     ylim([-50 5])
%     xlim([-FFTsize/2,FFTsize/2-1])
%     drawnow
%     legend('show')
%     legend('Location','Best')
% end

%% Clipping to improve PAPR
x_mag=abs(timeDomainOFDMSymbolsSerial);
x_max=0.6*max(x_mag);

for j=1:length(timeDomainOFDMSymbolsSerial)                                  %Clipping the signals above threshold(here 70% of original value)
if(x_mag(j)>x_max)
    x_clip(j)=x_max;
else
    x_clip(j)=x_mag(j);
end
end
Pmax = max((abs(x_clip)).^2);
Pavg = mean((abs(x_clip)).^2);

papr = Pmax/Pavg
plot(10*log10(abs(x_clip)))


%% Frame Generation and filtering
Pr = (mseq(2,6))';
Preambles = [Pr Pr];
% time_domain_symbols = time_domain_symbols/rms(time_domain_symbols);
Frame = [zeros(1,50) Preambles timeDomainOFDMSymbolsSerial zeros(1,50) ];
temp      = repmat(Frame,OSR,1);
TxFrame = reshape(temp,1,length(Frame)*OSR);


plot(20*log10((abs(TxFrame))))


% % Channel
[RxSignal,H_cir] = Channeling(TxFrame.',N_taps,fs);

waterfall(abs(H_cir))
plot(abs(RxSignal))

%% % % Synchronization
 [SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);
 scatterplot(SyncPreamble)

 estimated_channel = mean(SyncPreamble'./Preambles);

 receivedPreambles = SyncPreamble/estimated_channel;
 scatterplot(receivedPreambles)

 syncOFDMmatrix = reshape(SyncOFDMsignal,FFTsize+CPsize,[]);

 syncOFDMCPremoved = syncOFDMmatrix(CPsize+1:length(syncOFDMmatrix),:);

 fdomainSymbol = fft(syncOFDMCPremoved)/sqrt(FFTsize);

 receivedPilots = fdomainSymbol(pilotCarriers,:);

 estimatedOFDMChannel   = receivedPilots./pilotSymbols;
 
 receivedDataOFDM       = fdomainSymbol(availableCarriers,:);

% Channel Interpolation using interp
interPolated = interp1(pilotCarriers,estimatedOFDMChannel,1:1:length(fdomainSymbol),"spline");

interPolatedChannel = interPolated(availableCarriers,:);

equalized_frame = receivedDataOFDM./interPolatedChannel; 

equalized_frame(1) = [];
scatterplot(equalized_frame)
waterfall(abs(interPolatedChannel))
plot(abs(interPolated))


receivedBits = qamdemod(equalized_frame,2);
receivedBits(:)
receivedBits = [receivedBits zeros(1,9)];
%Binary to character conversion
receivedMessage = bin2char(receivedBits);




