% Get the serial number of the Pluto Radio
%PlutoRadioSerialNumbers = findPlutoRadio;

% Choose the first Pluto Radio as the transmitter
%TxRadioID = PlutoRadioSerialNumbers(1);
mod_stp.T='QAM';
mod_stp.M=4;
mod_stp.N=0.5e3;
flt_stp.T='RC';
flt_stp.sps=8;
flt_stp.span=12;
flt_stp.alpha=0.7;
flt_stp.BT=0.3;
flt_stp.mtchd= true;
hw_stp.Tx_G=-10;
hw_stp.Rx_G=20;
hw_stp.F=2.4e9;
hw_stp.R=1e6;
hw_stp.Tx_ID='sn:104473541196000cf7ff14008dd81e43f6';
hw_stp.Rx_ID='sn:1044735411960007f5ff2800134d1eb664';
hw_stp.N=30e3;

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
[frm, symbs_tx]=transmit(Tx,mod_stp,flt_stp);
[smpls_rx,prmb_rx,smpls_dt,symbs_rx]=receive(Rx,mod_stp,flt_stp,hw_stp);
release(Tx)
release(Rx)

% figure()
% 
% plot(10*log10(abs(smpls_rx).^2))
% grid on
% xlabel('time')
% ylabel('Power (dB)')
% title('power of the captured signal in time domain');
% figure()
% 
% plot(10*log10(abs(frm).^2))
% grid on
% xlabel('time')
% ylabel('Power (dB)')
% title('power of the transmitted signal in time domain');

figure()
[Pxx_rx,F_rx] = pwelch(smpls_rx,[],[],[],hw_stp.R,'centered');
plot(F_rx,10*log(Pxx_rx/max(Pxx_rx)));

% grid on
% xlabel('Frequency (Hz)')
% ylabel('Normalized Power (dB)')
% title('power spectrum of the captured signal');
% figure()
% [Pxx_tx,F_tx] = pwelch(frm,[],[],[],hw_stp.R,'centered');
% plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Normalized Power (dB)')
% title('power spectrum of the transmitted signal');
% hold on
% figure
%  obw(smpls_rx,hw_stp.R,[],99);
% 
% scatterplot(symbs_tx);
% 
% scatterplot(symbs_rx);

    eyediagram(symbs_rx16,4)
    % eyediagram(symbs_tx,6)
  %polar plots  
rx_signal  = symbs_rx
polarplot(angle(rx_signal),abs(rx_signal))
plot(real(rx_signal),imag(rx_signal));
grid on;
axis square;
xlabel('I-Component');
ylabel('Q-Component');
title('Rx signal');

ccdf = comm.CCDF('NumPoints', 10000000);
[ccdf_values, ccdf_bins] = step(ccdf, smpls_rx64);
figure()
semilogy(ccdf_bins, ccdf_values);
xlabel('Amplitude');
ylabel('Probability(%)');
title('Complementary Cumulative Distribution Function (CCDF)');
% figure()
% subplot(1,2,1)
% plot(real(smpls_rx),imag(smpls_rx));
% grid on;
% axis square;
% xlabel('I-Component');
% ylabel('Q-Component');
% title('Rx signal');
% 
% subplot(1,2,2)
% plot(real(frm),imag(frm));
% grid on;
% axis square;
% xlabel('I-Component');
% ylabel('Q-Component');
% title('Tx signal');
% 
% sgtitle('Polar Diagrams')
% evm = comm.EVM
% for i=1:length(symbs_rx)
% Evm(i) = evm(symbs_tx(i),symbs_rx(i))
% end
% Evm=Evm/length(symbs_rx);
% plot(Evm)
% ccdf = comm.CCDF
% 
% [ccx,ccy] = ccdf(smpls_rx);
% plot(ccx,ccy)