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
NumGuards            = 50;

% Message Source and Bits
Message = ['We hold these truths 2 be self-evident,' ...
    ' that all men are created equal, that they are ' ...
    'endowed by their Creator with certain unalienable ' ...
    'Rights, that among these are Life, Liberty and ' ...
    'the pursuit of Happiness'];

MessageBinChar      = dec2bin(double(char(Message)));
[Sz1,Sz2]           = size(MessageBinChar);
MessageBinStrn      = reshape(MessageBinChar,1,Sz1*Sz2);
MessageBits         = zeros(1,Sz1*Sz2);
for b=1:Sz1*Sz2
    MessageBits(b)  = str2double(MessageBinStrn(b));
end


% OFDM modulation
UsedSubcarrierIndx          = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];

pilotIndices                = UsedSubcarrierIndx(2:3:end);
UsedSubcarrierIndx(2:3:end) = [];

dataIndexes                 = UsedSubcarrierIndx;
pilotBits                   = repmat([1 0],1,length(pilotIndices)/2);
pilotSymbols                   = pilotBits*2-1;
%pilotSymbols                = pskmod(pilotBits,2);

numberOfDataCarriers        = NumActiveSubcarriers - length(pilotIndices); 
numberOfZerosToIncrease     = floor(length(MessageBits)/numberOfDataCarriers)*numberOfDataCarriers +numberOfDataCarriers - length(MessageBits);

ZeroPaddedMessageBits       = [MessageBits zeros(1,numberOfZerosToIncrease)];
modulatedDataSymbols        = pskmod(ZeroPaddedMessageBits,2);

ofdm_symbols_Matrix         = reshape(modulatedDataSymbols ,numberOfDataCarriers ,length(ZeroPaddedMessageBits)/numberOfDataCarriers);


to_ifft = zeros(FFTsize,length(ZeroPaddedMessageBits)/numberOfDataCarriers);

stem(abs(ofdm_symbols_Matrix))

for i= 1:length(ZeroPaddedMessageBits)/numberOfDataCarriers
    to_ifft(pilotIndices,i) = pilotSymbols;
end

to_ifft(dataIndexes,:)      = ofdm_symbols_Matrix;

ifftOFDM                    = ifft(to_ifft)*sqrt(FFTsize);

% iffOFDMtoSerial = reshape(ifftOFDM,1,length(ifftOFDM)*1);
ifftOFDMtoSerialwithCP      = [ifftOFDM((end-CPsize+1:end),:); ifftOFDM];
serialIFFT                  = ifftOFDMtoSerialwithCP(:).';

figure();
plot(10*log(abs(serialIFFT(1:FFTsize*4)).^2));
hold on
plot(10*log(abs(serialIFFT(1:FFTsize/CPsize)).^2));
hold off

figure()
[Pxx_tx_q_a,F_tx_q_a] = pwelch(serialIFFT,[],[],[],fs,'centered');
plot(F_tx_q_a,10*log10(Pxx_tx_q_a));

% mean max and PAPR

meanP = mean(abs(serialIFFT).^2);

peakP = max(abs(serialIFFT).^2);

PAPR2=peakP/meanP;

abs(serialIFFT).^2 >= meanP

% Frame Generation and filtering
preambles   =   mseq(2,6);
preambles = [preambles;preambles];

stem(mseq(2,6))
guardsymbols = zeros(1,NumGuards);

TxFrame = [guardsymbols preambles.' serialIFFT guardsymbols];
% TxFrame = [preambles;serialIFFT];

% temp      = repmat(TxFrame,OSR,1);
% TxFrame = reshape(temp,1,length(TxFrame)*OSR);

TxFrame = upsample(TxFrame,OSR);

plot(abs(TxFrame))

% Plotting TX Signal
figure()
plot(abs(TxFrame));

figure()
pwelch(TxFrame);

% % Channel
[RxSignal,H_cir] = Channeling(TxFrame.',N_taps,fs);
waterfall(abs(H_cir))

% Plotting Rx Signal
figure()
plot(abs(RxSignal));

figure()
plot(20*log10(abs(RxSignal)));


waterfall(abs(RxSignal))


% % Synchronization
 [SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);
 scatterplot(SyncPreamble)

% channel estimation
a = mean(SyncPreamble'./preambles.');
eq_pr = SyncPreamble/a;
scatterplot(eq_pr)

%% Windowing
%approach 1 Alphan's
w1 = 1/2*(cos(linspace(-pi,0,W+2))+1);
plot(w1)

w1(1) = [];
w1(end) = [];

W1 = eye(FFTsize+CPsize);
W1(1:W,1:W) = diag(w1);

w2 = 1/2*(cos(linspace(0,pi,W+2))+1);
plot(w2)

w2(1) = [];
w2(end) = [];
W2 = zeros(FFTsize+CPsize);

W2(1:W,(CPsize+(1:W))) = diag(w2);

xWOFDM = W1*ifftOFDMtoSerialwithCP  + W2* ifftOFDMtoSerialwithCP;

if (1) % OOB
    xA = ifftOFDMtoSerialwithCP(:);
    Nwin = 1024;
    Noverlap= floor(Nwin*3.2/4);
    [pX,f] = pwelch(xA(:), Nwin, Noverlap,[],FFTsize,'twosided');
    pX = pX * FFTsize;
    h100 = figure(100);
    plot(f-FFTsize/2,fftshift(10*log10(pX)) , 'r', 'displayname', 'OFDM' )
    hold on

    xA = xWOFDM(:).';
    Nwin = 1024;
    Noverlap= floor(Nwin*3.2/4);
    [pX,f] = pwelch(xA(:), Nwin, Noverlap,[],FFTsize,'twosided');
    pX = pX * FFTsize;
    h100 = figure(100);
    plot(f-FFTsize/2,fftshift(10*log10(pX)) , 'g', 'displayname', 'Windowing & OFDM' )
    hold on
    
    
    grid on
    xlabel('Subcarrier Index')
    ylabel('Power spectral density [dBm/Hz]')
    ylim([-50 5])
    xlim([-FFTsize/2,FFTsize/2-1])
    drawnow
    legend('show')
    legend('Location','Best')
end

%% Clipping 
Power=abs(TxFrame).^2
MeanPower=mean(abs(TxFrame).^2);

index=find(Power>=MeanPower)
tramsmitted_signalclip=TxFrame;

tramsmitted_signalclip(index)=10/100*tramsmitted_signalclip(index);

Papr=max(abs(TxFrame).^2)/mean(abs(TxFrame).^2)
Papr=max(abs(tramsmitted_signalclip).^2)/mean(abs(tramsmitted_signalclip).^2)


figure
plot(1:length(TxFrame),20*log(TxFrame))
hold on
plot(1:length(tramsmitted_signalclip),20*log(tramsmitted_signalclip))

figure
[Pxx_rx,F_rx] = pwelch(TxFrame,[],[],[],fs,'centered');
plot(F_rx,10*log10(Pxx_rx/max(Pxx_rx)))
hold on
[Pxx_rx,F_rx] = pwelch(tramsmitted_signalclip,[],[],[],fs,'centered');
plot(F_rx,10*log10(Pxx_rx/max(Pxx_rx)))

% OFDM demodulation, Channel estimation, equalization and detection

