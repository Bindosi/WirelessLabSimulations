clear; clc;

% Channel Parameters
Max_Excess_Delay   = 3.2*1e-5; % in seconds

% OFDM Signal Specs
FFTsize              = 64;
df                   = 15e3; % subcarrier spacing
fs                   = df*FFTsize;
Ts                   = 1/fs;
OSR                  = 1;
NumActiveSubcarriers = 48;
subCarriersPerUser   = 8;
numberOfUsers        = 4;
NumSymbols           = 10;
N_taps               = floor(Max_Excess_Delay/Ts)+1;
CPsizeOpt            = ceil(Max_Excess_Delay/Ts/OSR); % Minimum CP size required (effect of OSR taken into acount)
CPsize               = 8;
W                    = 2;
% loading the file
load("OFDMA_Data_setOSR1.mat");


UsedSubcarrierIndx          = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
availableCarriers           = UsedSubcarrierIndx;
pilotCarriers               = availableCarriers(2:3:end);
availableCarriers(2:3:end)  = [];
effectiveDataCarriers       = availableCarriers;

for userId = 1:numberOfUsers
    user(userId).Carriers           = effectiveDataCarriers(userId:3:end);
end
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
pilotBitMatrix         = repmat(pilotBits(:),1,NumSymbols);
pilotSymbols           = repmat(pilotBitsModulated(:),1,NumSymbols);


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
correlation         =   xcorr(PreamblesCreated,Y);
[~,indx1]           =   maxk(correlation,5);
sortedIndexes       =   sort(indx1);

figure()
plot(abs(correlation))

theActualFrameIndex = min(sortedIndexes)+28;
proposedframe       = Y(theActualFrameIndex :theActualFrameIndex+ ...
                     (CPsize+FFTsize)*NumSymbols +length(PreamblesCreated));

figure
plot(abs(proposedframe));
hold on
plot(abs(proposedframe(1:length(PreamblesCreated))))


figure()
[Pxx_tx,F_tx] = pwelch(proposedframe(127:end),[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))


prmbl_rx = proposedframe(1:length(PreamblesCreated)); 
phase    =  angle(sum(PreamblesCreated.*conj(prmbl_rx)));

equalize_preambles = prmbl_rx.*exp(j*phase);
scatterplot(equalize_preambles);


frameLength            = (FFTsize+CPsize)*NumSymbols;
symbols_rx             = proposedframe(length(prmbl_rx)+ ...
                         1:length(prmbl_rx)+frameLength);

%REMOVING CP
parallelOFDMSymbols    =   reshape(symbols_rx,FFTsize+CPsize,NumSymbols);
parOFDMSymbolsCPless   =   parallelOFDMSymbols(CPsize+1:end,:);

fdomainSymbol          = fft(parOFDMSymbolsCPless)/sqrt(FFTsize);
receivedPilots         = fdomainSymbol(pilotCarriers,:);
estimatedOFDMChannel   = receivedPilots./pilotSymbols;

% Channel Interpolation using interp
interPolated           = interp1(pilotCarriers,estimatedOFDMChannel,1:1:length(fdomainSymbol),"spline");
interPolatedChannel    = interPolated(availableCarriers,:);


receivedDataOFDM       = fdomainSymbol(availableCarriers,:);
scatterplot(receivedDataOFDM(:))

receivedDataOFDMFOC    = receivedDataOFDM(:)  .*exp(-j*phase);
scatterplot(receivedDataOFDMFOC)
equalized_frame        = receivedDataOFDM(:)./interPolatedChannel(:); 
scatterplot(equalized_frame)
xlim([-2 2])
ylim([-2 2])

receivedBits           = qamdemod(equalized_frame,2);
reshapeReceivedBits    = reshape(receivedBits,32,[]);
userBits               = zeros(FFTsize,NumSymbols);
userBits(effectiveDataCarriers,:) = reshapeReceivedBits;
receivedPilotBits      = qamdemod(receivedPilots,2);
% Getting bits for each user
for userId = 1:numberOfUsers
    user(userId).receivedBits     = userBits(user(userId).Carriers,:);
    user(userId).reshapedBits     = reshape(user(userId).receivedBits(1:77),7,[])';
    user(userId).messageBits      = user(userId).reshapedBits(:);
end

for userId                        = 1:numberOfUsers
    user(userId).Message          = char(reshape(bin2dec(reshape(char(user(userId).messageBits + '0'),7,[]).'),1,[]));
    msgbox(user(userId).Message,"The received Message")
end


BER_fading = length(find(receivedPilotBits(:) ~= pilotBitMatrix(:)))/length(pilotBitMatrix(:))



figure()
[Pxx_tx,F_tx] = pwelch(symbols_rx,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
xlabel('Frequency');
ylabel('Amplitude');
title('Frequency Spectrum of OFDM Symbols')
