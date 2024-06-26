% Written by Huseyin Arslan
% 1/05/2008
% Purpose: To teach students a basic digital communication simulation set up.
% Updated: 2/29/2024 by A. Kihero



%% DEFINE SIMULATION PARAMETERS
Rsym           = 1000;             % Symbol rate
OvSampRatio    = 16;               % Number of samples per symbol(Duration of the pulse)
Fs             = Rsym*OvSampRatio; % Sampling Frequency
ModOrder       = 4;                % QPSK
MessageLength  = 1000;             % Symbols
SNR            = 10;   

% Desired signal-to-noise ratio in [dB]
Values = [5];
for k =1:length(Values)
  SNR = Values(k);
    %% GENERATE MESSAGE BITS (random message)
NumBitsPerSym = log2(ModOrder); % Number of Bits Per Symbol
source_bits   = randn(1,MessageLength*NumBitsPerSym)>0; % Generate random bits-there are various other ways of doing the same thing



%%  MODULATION
if(0) % Our Own Implementation
    % Mapping the generated bits to QPSK symbols:
    % Angle [pi/4 3*pi/4 -3*pi/4 -pi/4] corresponds...
    % to Gray code vector [11 10 00 01], respectively.
    % 1 1 --> (+,+), 1 0 --> (+,-), 0 0 --> (-,-),  0 1 -->(-,+).

    table       = exp(1j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols
    table       = table([0 1 3 2]+1);      % Gray code mapping pattern for QPSK symbols
    inp         = reshape(source_bits,NumBitsPerSym,length(source_bits)/NumBitsPerSym);
    mod_symbols = table([2 1]*inp+1);      % maps transmitted bits into QPSK symbols
end


if(1) % MATLAB implementation using pskmod command
    inp         = reshape(source_bits,NumBitsPerSym,length(source_bits)/NumBitsPerSym);
    b2int       = [2 1]*inp; % this converts bits to integer index values (0,1,2,....7)
    mod_symbols = pskmod(b2int,ModOrder, pi/4);
end

if(0) % MATLAB's comm toolbox-based implementation
    myModulator = comm.PSKModulator(...
        'ModulationOrder',  ModOrder,...   % here the second input determines how to
        'PhaseOffset',      pi/4, ...      % rotate (phase offset) the constellation
        'BitInput',         true);         % point in the complex plane. If we don't
    % rotate, it provides QPSK costellations
    % on the I and Q axis like [1 j -1 -j].
    mod_symbols = transpose(myModulator(transpose(source_bits)));
end

% M = 16;
% OvSampRatio16QAM = 1;
% NumBitsPerSym16QAM          = log2(M); % Number of Bits Per Symbol
% source_bits_16QAM         = reshape(source_bits,NumBitsPerSym16QAM,length(source_bits)/NumBitsPerSym16QAM);
% dataSymbolsIn_16QAM = bit2int(source_bits_16QAM,NumBitsPerSym16QAM);
% mod_symbols_16QAM = qammod(dataSymbolsIn_16QAM,M); 
% EbNo = 10;
% snr = convertSNR(EbNo,'ebno', ...
%     samplespersymbol=OvSampRatio16QAM, ...
%     bitspersymbol=NumBitsPerSym16QAM);
% 
% receivedSignalG = awgn(mod_symbols_16QAM,snr,'measured');
% 
% scatterplot(mod_symbols_16QAM)
% grid on
% title('Constellation of the tx symbols 16QAM')
% 
% scatterplot(mod_symbols)
% grid on
% title('Constellation of the tx symbols')
% 
% sPlotFig = scatterplot(receivedSignalG,1,0,'g.');
% hold on
% scatterplot(mod_symbols_16QAM,1,0,'w*',sPlotFig)


%% PULSE SHAPING FILTERING
% for this experiment, we will use rectangular filter, later in the future
% labs, we will learn other types of filters

temp      = repmat(mod_symbols,OvSampRatio,1);
tx_signal = reshape(temp,1,length(mod_symbols)*OvSampRatio);

plotting real and imaginary parts for the transmitted signals
figure();
stem(real(tx_signal(1:6*OvSampRatio)))
hold on
plot(real(tx_signal(1:6*OvSampRatio)))
hold off


figure();
stem(imag(tx_signal(1:6*OvSampRatio)))
hold on
plot(imag(tx_signal(1:6*OvSampRatio)))
hold off
%% NOISE GENERATION FOR THE DESIRED SNR

len             = length(tx_signal);
noise_var = 10^(-SNR/10) ;
noise     = sqrt(noise_var) * (randn(1, len) + 1j*randn(1, len))/sqrt(2);



mean(real(noise))
mean(abs(noise))

plot(abs(noise))
figure;
hist3([real(noise);imag(noise)].','Nbins',[20 20])

figure;
plot(abs(xcorr(noise,noise))) 

%% PASSING THE MESSAGE THROUGH AWGN CHANNELsubplo

rx_signal = tx_signal + noise;

% calculation of actual signal and noise powers, and actual SNR value
signal_P = mean(abs(tx_signal).^2);
r_signal_P = mean(abs(rx_signal).^2);
noise_P  = mean(abs(noise).^2);
Current_true_SNR_Over_block_dB =10*log10(signal_P/noise_P);


figure;
subplot(1, 2, 1)
 histogram(real(noise),'Normalization','pdf')
axis square
title 'PDF of real'
grid on
subplot(1, 2, 2)
 histogram(imag(noise),'Normalization','pdf')
 title 'PDF of imaginary'
axis square
grid on

%Plotting imaginary parts of the received signal
figure;
plot(imag(tx_signal(1:6*OvSampRatio)))
axis square
title (sprintf('Imaginary Tx Signal SNR =  %d', SNR ))
grid on
hold on
plot(imag(rx_signal(1:6*OvSampRatio)))
 title ('Imaginary Received Signal')
axis square
grid on


%Plotting imaginary parts of the received signal
figure;
subplot(1, 2, 1)
plot(imag(tx_signal(1:6*OvSampRatio)))
axis square
title ('Imaginary Tx Signal')
grid on
subplot(1, 2, 2)
plot(imag(rx_signal(1:6*OvSampRatio)))
 title ('Imaginary Received Signal')
axis square
grid on

%Plotting imaginary parts of the received signal
figure;
subplot(1, 2, 1)
plot(real(tx_signal(1:6*OvSampRatio)))
axis square
title ('Real Tx Signal')
grid on
subplot(1, 2, 2)
plot(real(rx_signal(1:6*OvSampRatio)))
 title ('Real Received Signal')
axis square
grid on

% PLOTS
%% Power spectrum
figure()
[Pxx_tx,F_tx] = pwelch(tx_signal,[],[],[],Fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
grid on
xlabel('Frequency (Hz)')
ylabel('Normalized Power (dB)')
title('power spectrum of the Tx signal')

figure()
plot(rx_signal,'.')


%% Constellation diagram
if(1) % Our own Implementation
    figure()
    plot(rx_signal,'.')
    hold on
    plot(mod_symbols,'o','LineWidth',1.5)
    grid on;
    axis square;
    xlabel('I-Component');
    ylabel('Q-Component');
    title("Constellation of the Rx samples");
    legend('Rx', 'ideal');
    hold off
end

if(1) % Using Matlab scatterplot Command
    scatterplot(rx_signal)
    grid on
    title('Constellation of the Rx samples')
end
%% Polar diagram
figure()
subplot(1,2,1)
plot(real(tx_signal),imag(tx_signal));
grid on;
axis square;
xlabel('I-Component');
ylabel('Q-Component');
title('Tx signal');

subplot(1,2,2)
plot(real(rx_signal),imag(rx_signal));
grid on;
axis square;
xlabel('I-Component');
ylabel('Q-Component');
title('Rx signal');

sgtitle('Polar Diagrams')



%% Eye diagram
% Eye diagram
NumOfeyes = 8;
if(1) % Our own Implementation
    xTx = reshape(tx_signal,OvSampRatio*NumOfeyes,[]);
    figure
    subplot(2,1,1)
    plot(real(xTx(1:end-1,:)));
    grid on;
    xlabel('Time');
    ylabel('Amplitude');
    title("I-eye")

    subplot(2,1,2)
    plot(imag(xTx(1:end-1,:)));
    grid on;
    xlabel('Time');
    ylabel('Amplitude');
    title("Q-eye")
    sgtitle('Eye Diagrams: Tx signal')

    xRx = reshape(rx_signal,OvSampRatio*NumOfeyes,[]);
    figure
    subplot(2,1,1)
    plot(real(xRx(1:end-1,:)));
    grid on;
    xlabel('Time');
    ylabel('Amplitude');
    title("I-eye")

    subplot(2,1,2)
    plot(imag(xRx(1:end-1,:)));
    grid on;
    xlabel('Time');
    ylabel('Amplitude');
    title("Q-eye")

    sgtitle('Eye Diagrams: Rx signal')
end

if(0) % Using Matlab eyediagram Command
    eyediagram(tx_signal,4*OvSampRatio)
    eyediagram(rx_signal,4*OvSampRatio)
end



end



%% *********** BASIC RECEIVER IMPLEMENTATION
% This part will be developped by the students
