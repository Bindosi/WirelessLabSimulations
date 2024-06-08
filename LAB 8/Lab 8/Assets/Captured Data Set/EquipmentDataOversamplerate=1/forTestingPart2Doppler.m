clear; clc;

% Channel Parameters
Max_Excess_Delay   = 3.2*1e-8; % in seconds

% OFDM Signal Specs
FFTsize              = 64;
df                   = 5e3; % subcarrier spacing
fs                   = df*FFTsize;
Ts                   = 1/fs;
OSR                  = 1;
NumActiveSubcarriers = 48;
NumSymbols           = 25;
N_taps               = floor(Max_Excess_Delay/Ts)+1;
CPsizeOpt            = ceil(Max_Excess_Delay/Ts/OSR); % Minimum CP size required (effect of OSR taken into acount)
CPsize               = 4;
W                    = 2;

UsedSubcarrierIndx      = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
availableCarriers           = UsedSubcarrierIndx;
pilotCarriers               = availableCarriers(2:3:end);
availableCarriers(2:3:end)  = [];
effectiveDataCarriers   = availableCarriers;

% Generating Pilot Bits
pilotBits              = repmat([1 0],1, floor(length(pilotCarriers)/2));
% Modulating pilot bits and reshapinig them into pilotcarriers x Number of symbols
pilotBitsModulated     = pilotBits*2-1;
pilotSymbols           = repmat(pilotBitsModulated(:),1,NumSymbols);

% loading the file
load("OFDM_RVC_CPSize4.mat");
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
correlation         =   xcorr(Y,PreamblesCreated)
[~,indx1]           =   maxk(correlation,5)

sortedIndex         = sort(indx1)
frameStartIndex     = sortedIndex(2) - length(Y) -1;
proposedframe = Y(frameStartIndex:frameStartIndex+length(PreamblesCreated)+ (CPsize+FFTsize)*NumSymbols)

figure()
plot(abs(correlation))

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

receivedDataOFDM        = fdomainSymbol(availableCarriers,:);
equalized_frame         = receivedDataOFDM(:)./interPolatedChannel(:); 
receivedBits            = qamdemod(equalized_frame,2);

scatterplot(equalized_frame); 


figure()
[Pxx_tx,F_tx] = pwelch(symbols_rx,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))

symbols_off=1
    indx = 1;
    prmbl_length = length(PreamblesCreated);
    crr_prv = 0;
    for i = 1:2*symbols_off
        crr = abs(PreamblesCreated.'*proposedframe(i:i+prmbl_length-1));
        correl(i) = crr; 
        if crr > crr_prv
            crr_prv = crr;
            indx = i;
        end
        
    end
    figure()
    plot(correl);
    title('Correlation for the preambles')
    indx_start  = indx+prmbl_length;




