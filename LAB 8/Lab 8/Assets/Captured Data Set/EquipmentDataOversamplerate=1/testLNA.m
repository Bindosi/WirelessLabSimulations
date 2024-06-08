clear; clc;

% Channel Parameters
Max_Excess_Delay   = 3.2*1e-5; % in seconds

% OFDM Signal Specs
FFTsize              = 64;
df                   = 15e3; % subcarrier spacing
fs                   = df*FFTsize;
sys.fs               = fs
Ts                   = 1/fs;
OSR                  = 1;
sys.sps              = OSR;
NumActiveSubcarriers = 32;
NumSymbols           = 5;
N_taps               = floor(Max_Excess_Delay/Ts)+1;
CPsizeOpt            = ceil(Max_Excess_Delay/Ts/OSR); % Minimum CP size required (effect of OSR taken into acount)
CPsize               = 8;
W                    = 2;

UsedSubcarrierIndx          = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
availableCarriers           = UsedSubcarrierIndx;
pilotCarriers               = availableCarriers(2:3:end);
availableCarriers(2:3:end)  = [];
effectiveDataCarriers   = availableCarriers;

% Generating Pilot Bits
pilotBits              = repmat([1 0],1,floor(length(pilotCarriers)/2));

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

% loading the file
load("OFDM_LNA_S32.mat");
% plot(abs(Y))
figure
plot(abs(Y))

Freq=fftshift(fft(Y))

figure
plot(abs(Freq))
 figure()
[Pxx_tx,F_tx] = pwelch(Y,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))

%
PreamblesCreated    =   [mseq(2,6); mseq(2,6)];
correlation         =   xcorr(PreamblesCreated,Y)
[~,indx1]           =   maxk(correlation,5)

figure()
plot(abs(correlation))

theActualFrameIndex      = length(Y) - indx1
proposedframe = Y(theActualFrameIndex :theActualFrameIndex+ (CPsize+FFTsize)*NumSymbols +length(PreamblesCreated))

figure
subplot(2,1,1),plot(abs(proposedframe))
subplot(2,1,2),plot(1:length(proposedframe),20*log10(proposedframe))

figure()
[Pxx_tx,F_tx] = pwelch(proposedframe(127:end),[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))


prmbl_rx = proposedframe(1:length(PreamblesCreated)); 
phase    =  angle(sum(PreamblesCreated.*conj(prmbl_rx)));

equalize_preambles = prmbl_rx.*exp(j*phase);
scatterplot(equalize_preambles);

frameLength = (FFTsize+CPsize)*NumSymbols;
symbols_rx = proposedframe(length(prmbl_rx)+1:length(prmbl_rx)+frameLength);

%REMOVING CP
parallelOFDMSymbols     =   reshape(symbols_rx,FFTsize+CPsize,NumSymbols)
parOFDMSymbolsCPless    =   parallelOFDMSymbols(CPsize+1:end,:)

fdomainSymbol          = fft(parOFDMSymbolsCPless)/sqrt(FFTsize);
receivedPilots         = fdomainSymbol(pilotCarriers,:);
estimatedOFDMChannel   = receivedPilots./pilotSymbols;

% Channel Interpolation using interp
interPolated            = interp1(pilotCarriers,estimatedOFDMChannel,1:1:length(fdomainSymbol),"spline");
interPolatedChannel     = interPolated(availableCarriers,:);


receivedDataOFDM           = fdomainSymbol(availableCarriers,:);
scatterplot(receivedDataOFDM(:))

if(1)
    symbls_off = 10;
    coarse_offset = est_coarse_fre(prmbl_rx, sys);
    % Compensation of frequency offset
    t = (1:length(receivedDataOFDM(:)))'/sys.fs;
    % compensation of the frequency offset
    receivedDataOFDM = receivedDataOFDM(:).*exp(-2i*pi*t*coarse_offset);
  
          scatterplot(receivedDataOFDM)
         % eyediagram(smpls_rx_q,256)
end
receivedDataOFDMFOC        = receivedDataOFDM;
% scatterplot(receivedDataOFDMFOC)
equalized_frame         = receivedDataOFDMFOC./interPolatedChannel(:); 
receivedBits            = qamdemod(equalized_frame,2);

scatterplot(equalized_frame);
xlim([-2 2])
ylim([-2 2])


figure()
[Pxx_tx,F_tx] = pwelch(symbols_rx,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))


function off = est_coarse_fre(preambles,sys)
    
    sampx = fftshift(fft(preambles.^2,40*length(preambles)));
    delta_f = (sys.fs)/length(sampx);

    freq_interv=(-(sys.fs))/2:delta_f:((sys.fs)/2)-1;
    figure()
    plot(freq_interv/2, abs(sampx));
    off  = freq_interv(sampx==max(sampx))/2;
end