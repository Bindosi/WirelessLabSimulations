mod_stp.T = 'QAM';
mod_stp.M = 16;
mod_stp.N = 0.5e3;


flt_stp.T = 'RC';
flt_stp.sps = 8;
flt_stp.span = 12;
flt_stp.alpha = 0.3;
flt_stp.BT = 0.7;
flt_stp.mtchd = true;

%Hardware Setup
hw_stp.Tx_G  = -20;
hw_stp.Rx_G  = 20;
hw_stp.F     = 2.4e9;
hw_stp.R     = 1e6;
hw_stp.Tx_ID = 'sn:10447354119600060d001800cf281e583b';
hw_stp.Rx_ID ='sn:104473dc599300131200210082672a4170';
hw_stp.N = 40e3;

Tx = sdrtx('Pluto', 'Gain',hw_stp.Tx_G,...
    'CenterFrequency',hw_stp.F,...
    'BasebandSampleRate',hw_stp.R,...
    'RadioID',hw_stp.Tx_ID);

Rx = sdrrx('Pluto', 'SamplesPerFrame', hw_stp.N,...
    'OutputDataType','double',...
    'CenterFrequency',hw_stp.F,...
    'BasebandSampleRate',hw_stp.R,...
    'GainSource','Manual',... %AGC Fast Attack
    'Gain', hw_stp.Rx_G,...
    'RadioID',hw_stp.Rx_ID);

%Transmiting and receiving
[frm,symbs_tx] = transmit(Tx, mod_stp, flt_stp);
[smpls_rx,prmb_rx,smpls_dt,symbs_rx] = receive(Rx,mod_stp, flt_stp,hw_stp);

% smpls_rx = received raw data samples,
% prmb_rx  = received preamble samples,
% smpls_dt  = received processed data samples,
% symbs_rx = received data samples

[Pxx_tx,F_tx] = pwelch(smpls_rx,[],[],[],hw_stp.R,'centered');


% figure()
% plot(20*log10(abs(smpls_rx)));
% grid on
% xlabel('Time');
% ylabel('Normalized Power (dB)')
% title('Power of received signal in time domain')

% % polar diagram
% plot(real(smpls_dt),imag(smpls_dt));
% grid on;
% axis square;
% xlabel('I-Component');
% ylabel('Q-Component');
% title('Rx signal');
% 
% sgtitle('Polar Diagrams 0.3 alpha (roll off factor)')

%99 and 50% ocuupied bandwidth
% figure
obw(smpls_rx,hw_stp.R,[],99);

ccdf = comm.CCDF(...
    'AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true);
[ccdfy,ccdfx,avg,peak] = ccdf([smpls_rx smpls_rx]);

pm = powermeter(ComputeCCDF=true);
averagePower = pm(smpls_rx); 
prob = probability(pm,3)
plotCCDF(pm,GaussianReference=true)
 title('CCDF Measurements 16QAM')

plot(ccdf)
legend('16-QAM','QPSK')



%constellation diagrams
    scatterplot(symbs_rx)
    grid on
    title('Constellation of the Rx samples')
%EVM machine
    plot(comm.EVM(symbs_tx,smpls_rx ))

tscope = timescope(YLabel="EVM (%)",YLimits=[0 40], ...
    SampleRate=1000,TimeSpanSource="Property",TimeSpan=1.2, ...
    ShowGrid=true);
evm = comm.EVM;
for i = 1:length(symbs_tx)
e(i) = evm(symbs_tx(i),symbs_rx(i));
end
e=e/length(symbs_tx)
figure
plot(e)
title('Error Vector Magnitude for gain = -10')
xlabel('Time')
ylabel('EVM')
tscope(e)
% %Eye diagrams
% eyediagram(smpls_dt, 16)
% eyediagram(frm, 16)

% plot(F_tx,10*log10(Pxx_tx/max(Pxx_tx)))
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Normalized Power (dB)')
% title('Power spectrum of the Rx signal Part II  Q4')

release(Tx);
release(Rx);