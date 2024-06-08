% Wireless Lab Course
% Istanbul Medipol University
% CoSiNC Research Group
% Spring 2023/2024

% Lab #1: SIMULATION OF A SIMPLE COMMUNICATION LINK WITH FREQUENCY OFFSET.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
% SIMULATION PARAMETERS
Rsym           = 10;                % Symbol rate
OvSampRatio    = 256;               % Number of samples per symbol(Duration of the pulse)
Fs             = Rsym*OvSampRatio;  % Sampling Frequency
ModOrder       = 4;                 % QPSK
MessageLength  = 200;               % Bits
SNR            = 100;               % Desired signal-to-noise ratio in [dB]
F_off          = [25 0 -115];       % Frequency Offset Values in Hz.

% GENERATING RANDOM MESSAGE BITS
% there are various other ways of doing the same thing
source_bits  = randn(1,MessageLength) > 0;  % Message bits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QPSK MODULATION
% Mapping the generated bits to QPSK symbols:
% Angle [pi/4 3*pi/4 -3*pi/4 -pi/4] corresponds...
% to Gray code vector [11 10 00 01], respectively.
% 1 1 --> (+,+), 1 0 --> (+,-), 0 0 --> (-,-),  0 1 -->(-,+).

if(1)
    table       = exp(1j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols
    table       = table([0 1 3 2]+1);      % Gray code mapping pattern for QPSK symbols
    inp         = reshape(source_bits,2,length(source_bits)/2);
    mod_symbols = table([2 1]*inp+1);      % maps transmitted bits into QPSK symbols
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PULSE SHAPING FILTERING
% for simplicity, we use Rectangular filter here
temp      = repmat(mod_symbols,OvSampRatio,1);
tx_signal = reshape(temp,1,length(mod_symbols)*OvSampRatio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AWGN CHANNEL
% Noise generation for the desired SNR value
len       = length(tx_signal);
noise_var = 10^(-SNR/10);
noise     = sqrt(noise_var)*(randn(1, len) + 1j*randn(1, len))/sqrt(2);

% Rx Signal
rx_signal = tx_signal + noise;

% Introduce Frequency offset.
sampleIndices = 0:length(rx_signal)-1; % in samples
t             = sampleIndices./Fs;     % in second
for ii=1:length(F_off)
    Phase_offset(ii,:)   = 2*pi*F_off(ii)*t;
    CarrierF_off(ii,:)   = exp(1j*Phase_offset(ii,:));
    rx_signal_Foff(ii,:) = rx_signal.*CarrierF_off(ii,:);
end

%%%%%%%%%%%%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First 10 symbols
figure
for ii = 1:length(F_off)
    temp = F_off(ii);
    subplot(2,3,(ii))
    plot(real(rx_signal_Foff(ii,1:6*OvSampRatio))); axis square
    xlabel('Time'); ylabel('Amplitude'); title("Real part: Foff =" + temp+ "Hz")
    subplot(2,3,(ii+3))
    plot(imag(rx_signal_Foff(ii,1:6*OvSampRatio))); axis square
    xlabel('Time'); ylabel('Amplitude');title("Imaginary part: Foff =" + temp+ "Hz")
end
sgtitle('Time Domain Rx signals')

% Power spectrum
legendCell = cellstr(num2str(F_off', 'Foff=%-d'));
legendCell = ['TxSignal'; legendCell];

figure
[Pxx_tx,F_tx] = pwelch(tx_signal,[],[],[],Fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
grid on; box on
xlabel('Frequency (Hz)')
ylabel('Normalized Power (dB)')
title('power spectrum of the Tx and Rx signals')
hold on
for ii = 1:length(F_off)
    [Pxx_rx,F_rx] = pwelch(rx_signal_Foff(ii,:),[],[],[],Fs,'centered');
    plot(F_rx,10*log10(Pxx_rx/max(Pxx_rx)))
end
legend(legendCell)
hold off

% Constellation diagram
figure
for ii = 1:length(F_off)
    temp = F_off(ii);
    subplot(1,3,ii)
    plot(rx_signal_Foff(ii,:),'.')
    hold on
    plot(mod_symbols,'o','LineWidth',1.5)
    grid on; axis square
    title("constellation: Foff=" + temp+"Hz")
    legend('Rx', 'ideal')
    hold off
end
sgtitle('Constellation Diagrams')

% Polar diagram
figure
for ii = 1:length(F_off)
    temp = F_off(ii);
    subplot(1,3,ii)
    plot(real(rx_signal_Foff(ii,:)),imag(rx_signal_Foff(ii,:)));
    grid on; axis square
    title("polar Diagram: Foff=" + temp+"Hz")
end
sgtitle('Polar Diagrams')

% Eye diagram
figure
for ii = 1:length(F_off)
    x   = reshape(rx_signal_Foff(ii,:),OvSampRatio*4,[]);
    temp = F_off(ii);
    subplot(2,3,(ii))
    plot(real(x(1:end-1,1:10))); axis square
    xlabel('Time'); ylabel('Amplitude'); title("I-eye: Foff =" + temp+ "Hz")
    subplot(2,3,(ii+3))
    plot(imag(x(1:end-1,1:10))); axis square
    xlabel('Time'); ylabel('Amplitude');title("Q-eye: Foff =" + temp+ "Hz")
end
sgtitle('Eye Diagrams')


