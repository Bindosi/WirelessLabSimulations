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

% Message Source and Bits
Message = ['We hold these truths 2 be self-evident,' ...
    ' that all men are created equal, that they are ' ...
    'endowed by their Creator with certain unalienable ' ...
    'Rights, that among these are Life, Liberty and ' ...
    'the pursuit of Happiness'];

MessageBinChar = dec2bin(double(char(Message)));
[Sz1,Sz2]      = size(MessageBinChar);
MessageBinStrn = reshape(MessageBinChar',1,Sz1*Sz2);
MessageBits    = zeros(1,Sz1*Sz2);
for b=1:Sz1*Sz2
    MessageBits(b) = str2double(MessageBinStrn(b));
end


% % 
%Pilot bit generation and insertiom
pilots      = [0 1];
FrameBits   = [];
k = 0;
for i=1:1+(length(MessageBits)-1)/2+length(MessageBits)
k=k+1;
    if rem(i,3)==2
           ind = ((-1)^(i)+1)/2+1;
    FrameBits(i) = pilots(ind);
     k=k-1;
    else
    FrameBits(i)=MessageBits(k);

    end
end


% temp = [0 MessageBits];
% q = length(temp)/2;
% temp = reshape(temp,q,2);
% pilots = zeros(q,1);
% pilots(1:2:end-1)=1;
% temp = [temp pilots];
% temp= reshape(temp',1,[]);
% FrameBits =temp(2:end);


%BPSK
symbols = FrameBits*2-1;

% OFDM modulation


NumActiveSubcarriers = 48;
UsedSubcarrierIndx = [FFTsize-(NumActiveSubcarriers/2:-1:1)+1, 1:NumActiveSubcarriers/2];
temp1 = rem(length(symbols),NumActiveSubcarriers);

symbols = [symbols zeros(1,temp1)]; % 14 being the remainder of length(symbols) over 48;
time_domain_symbols = [];
nSyms = 0;
for i=1:NumActiveSubcarriers:(length(symbols)-NumActiveSubcarriers+1)
nSyms=nSyms+1;
    temps = symbols(i:i+NumActiveSubcarriers-1);
    freq = zeros(1,64);
    freq(UsedSubcarrierIndx)=temps;
      time = ifft(freq);%*sqrt(FFTsize);
    time = [time(end-CPsize+1:end) time];
time_domain_symbols = [time_domain_symbols time];
end

time_domain_symbols = time_domain_symbols(1:25*(64+CPsize));
length(FrameBits)
figure()
plot((abs(time_domain_symbols)))
% PAPR 
Pmax = max((abs(time_domain_symbols)).^2);
Pavg = mean((abs(time_domain_symbols)).^2);

papr = Pmax/Pavg

%CCDF
ccdf = comm.CCDF('NumPoints', 1000000);
[ccdf_values, ccdf_bins] = step(ccdf, time_domain_symbols');
figure()
semilogy(ccdf_bins, ccdf_values);
xlabel('Amplitude');
ylabel('Probability(%)');
title('Complementary Cumulative Distribution Function (CCDF)');
%% Enhancing OFDM


%spectrum assuming a symbol rate of fs

figure()
[Pxx_tx,F_tx] = pwelch(time_domain_symbols,[],[],[],fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
%% Clipping
for i =1:length(time_domain_symbols)
    if abs(time_domain_symbols(i))^2/Pavg > .5*papr
        clip(i) = sqrt(Pavg)*time_domain_symbols(i)/abs(time_domain_symbols(i));
    else
        clip(i) = time_domain_symbols(i);
    end
end

Pmax = max((abs(clip)).^2);
Pavg = mean((abs(clip)).^2);

papr = Pmax/Pavg

plot(10*log10(abs(clip)))


%% Frame Generation and filtering
Pr = (mseq(2,6))';
% time_domain_symbols = time_domain_symbols/rms(time_domain_symbols);
Frame = [zeros(1,50) Pr Pr time_domain_symbols zeros(1,50)];
temp      = repmat(Frame,OSR,1);
TxFrame = reshape(temp,1,length(Frame)*OSR);


plot(20*log10((abs(TxFrame))))


% % Channel
[RxSignal,H_cir] = Channeling(TxFrame.',N_taps,fs);

waterfall(abs(H_cir))
plot(abs(RxSignal))
length(time_domain_symbols)*OSR

% % Synchronization

[SyncPreamble,SyncOFDMsignal] = synchronization(RxSignal,OSR,fs,FFTsize,CPsize,NumSymbols);

scatterplot(SyncPreamble)

% channel estimation
a = mean(SyncPreamble'./[Pr Pr]);
eq_pr = SyncPreamble/a;
scatterplot(eq_pr)
% h = zeros(1,N_taps);
% h(1) = SyncPreamble(1)/Pr(1);
% for i = 1:31
%    % h(1)*Pr(i)+h(2)*Pr(i-1)+...+h(i)*Pr(1)=SyncPreamble(i)
%         h(i)=(SyncPreamble(i)-sum((Pr(i:-1:2)).*(h(1:1:i-1))))/Pr(1);
% end

%% OFDM demodulation, Channel estimation, equalization and detection
ofdm = [];
jmp = 64+CPsize;

for i=1:jmp:length(SyncOFDMsignal)-jmp+1
    to_fft= SyncOFDMsignal(i+CPsize:i+64+CPsize-1);
   temp = fft(to_fft);
   ofdm = [ofdm temp'];
end

length(temp)
scatterplot(ofdm)

% using pilots to estimate

h1 = ofdm(2:3:1598)./symbols(2:3:1598);
% we'll have 1 symbol at the beginning and two at the end
vq2 = interp1(2:3:1598,h1,2:1:1598,'spline');
h= [1 vq2 1 1];
eq = ofdm(2:1598)./vq2;
scatterplot(eq/rms(eq));

det = (sign(eq)+1)/2;
bits = [];
bits = [det(1) bits];
k = 1;
for i=3:3:length(det)
    k= k+1;
bits(k) = det(i);
bits(k+1) = det(i+1);

