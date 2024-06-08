Fs = 102.4e3;
% loading signals from dataset
load ("sig_hops.mat")   % loading signal wave
Set1Signal1_Y       = loadStructFromFile("Set1Sig1.mat", "Y"); % loading Y into Set1Signal1_Y  
Set1Signal2_Y       = loadStructFromFile("Set1Sig2.mat", "Y"); % loading Y into Set1Signal2_Y 
pspectrum(wave)
load("Set2Sig1.mat")

Set2Sig1            = loadStructFromFile("Set2Sig1.mat", "FilteredFrame");  
Set2Sig2            = loadStructFromFile("Set2Sig2.mat", "FilteredFrame");  
Set2Sig3            = loadStructFromFile("Set2Sig3.mat", "FilteredFrame");  

% Plotting the dynamic ranges of signals
subplot(3,1,1)
envelope(abs(Set2Sig1(100:7500)))
title('Analytic Envelopg Signal Set2Sig1')
% [Set2Sig1_up, Set2Sig1_lo] = envelope(abs(Set2Sig1));
% hold on
% plot(Set2Sig1,Set2Sig1_up, Set2Sig1, Set2Sig1_lo,'linewidth',1.5)
% hold off
subplot(3,1,2)
envelope(abs(Set2Sig2(100:7500)))
title('Analytic Envelopg Signal Set2Sig2')
% [Set2Sig2_up,Set2Sig2_lo] = envelope(abs(Set2Sig2));
% hold on
% plot(Set2Sig2,Set2Sig2_up, Set2Sig2, Set2Sig2_lo,'linewidth',1.5)
% hold off
subplot(3,1,3)
envelope(abs(Set2Sig3(100:7500)))
title('Analytic Envelopg Signal Set2Sig3')
% [Set2Sig3_up,Set2Sig3_lo] = envelope(abs(Set2Sig3));
% hold on
% plot(Set2Sig3, Set2Sig3_up, Set2Sig3, Set2Sig3_up,'linewidth',1.5)
% hold off


% Plotting the polar plots of signals

subplot(3,1,1)
plot(real(Set2Sig1),imag(Set2Sig1));
grid on;
axis square;
xlabel('I-Component');
ylabel('Q-Component');
title('Set2Sig1 Polar Diagram');

subplot(3,1,2)
plot(real(Set2Sig2),imag(Set2Sig2));
grid on;
xlabel('I-Component');
ylabel('Q-Component');
title('Set2Sig2 Polar Diagram');

subplot(3,1,3)
plot(real(Set2Sig3),imag(Set2Sig3));
grid on;
axis square;
xlabel('I-Component');
ylabel('Q-Component');
title('Set2Sig3 Polar Diagram');

figure()
plot(10*log10(abs(smpls_rx)));
xlabel('Time')
ylabel('Normalized Power (dB)')
title('power spectrum of time domain of the Rx signal Wireless with QAM')

% plotting for the first signal
figure()
[Pxx_tx,F_tx] = pwelch(Set1Signal1_Y,[],[],[], Fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
grid on
xlabel('Frequency (Hz)')
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Set1Signal1 Signal')

% plotting for the second signal
figure()
[Pxx_tx,F_tx]   = pwelch(Set1Signal2_Y,[],[],[], Fs,'centered');
plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
grid on
xlabel('Frequency (Hz)')
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Set1Signal2 Signal')

%Plotting CCDF curves for matlab signal
pm_1            = powermeter(ComputeCCDF=true);
pm_2            = powermeter(ComputeCCDF=true);
pm_3            = powermeter(ComputeCCDF=true);

figure()
averagePower    = pm_1(Set2Sig1);
plotCCDF(pm_1)
hold on

averagePower    = pm_2(Set2Sig2);
plotCCDF(pm_2)
hold on

averagePower    = pm_3(Set2Sig3);
plotCCDF(pm_3)
legend('Sig2Sig1','Sig2Sig2','Sig2Sig3')
hold off

% Plotting eye diagrams

eyediagram(Set2Sig1,16)


eyediagram(Set2Sig2,16)


eyediagram(Set2Sig3,16)


