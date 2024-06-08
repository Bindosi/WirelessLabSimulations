
% Filter Setup
s.T ="RC";
s.sps = 8;
s.span = 12;
s.alpha = 0.5;
%Preamble
Pr1 = mseq(2,7);
Pr2 = mseq(2,7);
%Symbols PSK
source_bits = (randn(1,2*512))>0;
myModulator = comm.PSKModulator(...
        'ModulationOrder',  4,...  
        'PhaseOffset',      pi/4, ...     
        'BitInput',         true);         
 mod_symbols = transpose(myModulator(transpose(source_bits)));
hw_stp.Tx_G=0;
hw_stp.Rx_G=20;
hw_stp.F=2.4e9;
hw_stp.R=1e6;
hw_stp.Tx_ID='sn:10447354119600060d001800cf281e583b';
hw_stp.Rx_ID='sn:104473dc599300131200210082672a4170';
hw_stp.N=40e3;

% Initialize Tx object with the given parameters
Tx = sdrtx( 'Pluto','Gain', hw_stp.Tx_G, ...
    'CenterFrequency', hw_stp.F, ...
    'BasebandSampleRate',hw_stp.R, ...
    'RadioID', hw_stp.Tx_ID);
% Choose the second Pluto Radio as the receiver
%RxRadioID = PlutoRadioSerialNumbers(2);

% Initialize Rx object with the given parameters
Rx = sdrrx( 'Pluto', 'SamplesPerFrame', hw_stp.N, ...
    'OutputDataType', 'double', ...
    'CenterFrequency', hw_stp.F, ...
    'BasebandSampleRate',hw_stp.R,...
    'GainSource', 'Manual', ...
    'Gain', hw_stp.Rx_G, ...
    'RadioID', hw_stp.Rx_ID);
temp = [Pr1; Pr2; mod_symbols'];
frm1Tx = fltr(temp, s);
mod_stp.T='QPSK';
mod_stp.M=4;

% NOT UNDERSTOOD WHY MOD:STP.N = 512
mod_stp.N=512;
Tx.transmitRepeat(frm1Tx); 
[smpls_rx1,prmb_rx1,smpls_dt1,symbs_rx1]=receive(Rx,mod_stp,s,hw_stp);



release(Tx)
release(Rx)



figure()
[Pxx_tx,F_tx] = pwelch(smpls_rx1,[],[],[],hw_stp.R,'centered');
plot(F_tx,10*log10(Pxx_tx));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Rx signal')

figure()
obw(smpls_rx1,hw_stp.R,[],99);
P_rx1 = mean(Pxx_tx);

% Question 6
e = comm.EVM
e = comm.EVM
evm2 = e(temp(end-511:end),symbs_rx1/rms(symbs_rx1));

evm2 = [7.13 4.92 7.22 7.35 5.49 7.37]
SNR = -20*log10(evm2/100);
SNR = SNR'



%question 7 noise

[smpls_rx7,prmb_rx1,smpls_dt1,symbs_rx1]=receive(Rx,mod_stp,s,hw_stp);
release(Rx)
figure()
x = xcorr(smpls_rx7-mean((smpls_rx7)));
Xf = fftshift(fft(x));
% [Pxx_tx,F_tx] = pwelch(smpls_rx1,[],[],[],hw_stp.R,'centered');
[Xf,F_tx] = pwelch(smpls_rx7,[],[],[],hw_stp.R,'centered');
plot(F_tx,10*log10(Xf));
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Noise signal')
P_rx7 = mean(Xf);

totPower2 = trapz(Xf);
plot(10*log((abs(x)).^2))

%% Solving the Questions

%5. Plot the power spectrum of the received signal and use it to calculate the received signal power.

figure()
[Pxx_tx,F_tx] = pwelch(smpls_rx1,[],[],[],hw_stp.R,'centered');
plot(F_tx,10*log10(Pxx_tx));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)');
title('power spectrum of the Rx signal');
obw(smpls_rx1,hw_stp.R,[],99);
P_rx1 = mean(Pxx_tx);
[Pxx, F] = pspectrum(smpls_rx1, hw_stp.R);
plot(F,Pxx);
totPower = trapz(abs(Pxx_tx), F_tx);
totPower2 = trapz(abs(smpls_rx1).^2);

% we still need to calculate the power and the EVM
% from Q6
e = comm.EVM
evm1temp = e(symbs_rx1,OriginalSyms);

% scatterplot(mod_symbols_temp);
% 
% scatterplot(temp(end-511:end));
% mod_symbols_temp = [];
% 
% for ii=1:length(symbs_rx1)
%     if (real(symbs_rx1(ii))>0)
%         if (imag(symbs_rx1(ii))>0)
%             mod_symbols_temp(ii) = exp(j*pi*1/4);
%         else 
%             mod_symbols_temp(ii) = exp(-j*pi*1/4);
%         end 
%     else
%         if (imag(symbs_rx1(ii))>0)
%             mod_symbols_temp(ii) = exp(j*pi*3/4);
%         else 
%             mod_symbols_temp(ii) = exp(-j*pi*3/4);
%         end 
%     end
% end
% mod_symbols_temp = zeros(size(symbs_rx1)); 
% 
% for ii = 1:length(symbs_rx1)
%     if (real(symbs_rx1(ii)) > 0)
%         if (imag(symbs_rx1(ii)) > 0)
%             mod_symbols_temp(ii) = exp(j*pi*1/4);
%         else 
%             mod_symbols_temp(ii) = exp(-j*pi*1/4);
%         end 
%     else
%         if (imag(symbs_rx1(ii)) > 0)
%             mod_symbols_temp(ii) = exp(j*pi*3/4);
%         else 
%             mod_symbols_temp(ii) = exp(-j*pi*3/4);
%         end 
%     end
% end
% mod_symbols_temp = mod_symbols_temp';
% mod_symbols_temp9 = [];
% 
% for ii=1:length(symbs_rx91)
%     if real(symbs_rx91(ii))>0
%         if imag(symbs_rx91(ii))>0
%             mod_symbols_temp9(ii) = exp(j*pi*1/4);
%         else 
%             mod_symbols_temp9(ii) = exp(-j*pi*1/4);
%         end 
%     else
%         if imag(symbs_rx91(ii))>0
%             mod_symbols_temp9(ii) = exp(j*pi*3/4);
%         else 
%             mod_symbols_temp9(ii) = exp(-j*pi*3/4);
%         end 
%     end
% end
% scatterplot(symbs_rx1/rms(symbs_rx1));
% scatterplot(mod_symbols_temp);
% scatterplot(temp(end-511:end));
% scatterplot(symbs_rx91/rms(symbs_rx91));
%7  the noise
figure()
x = xcorr(smpls_rx7-mean((smpls_rx7)));
Xf = fftshift(fft(x));
plot(10*log((abs(x)).^2));
P_rx7 = mean(Pxx_tx);

% Question 9 increasing the Gain 

hw_stp.Tx_G=0;

e = comm.EVM
evm2 = e(temp(end-511:end),symbs_rx91/rms(symbs_rx91));
%% Playing Around

% observing what happens to the sig at every stage
% befor filter: mod_Symbols
scatterplot(mod_symbols);
plot(abs(mod_symbols));
% After adding preamble temp
scatterplot(temp);
plot(abs(temp));
stem(linspace(1,40000,length(temp)),abs(temp));
% after the filter frm1Tx
scatterplot(frm1Tx);
stem(abs(frm1Tx));
