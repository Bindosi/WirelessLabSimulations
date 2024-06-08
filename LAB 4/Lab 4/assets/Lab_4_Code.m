mod_stp.T = 'QPSK';
mod_stp.M = 4;
mod_stp.N = 512;


% flt_stp.T = 'RRC';
flt_stp.sps = 8;
flt_stp.span = 12;
flt_stp.alpha = 0.3;
% flt_stp.BT = 0.7;
flt_stp.mtchd = true;

%Hardware Setup
hw_stp.Tx_G  = -10;
hw_stp.Rx_G  = 20;
hw_stp.F     = 2.4e9;
hw_stp.R     = 1e6;
hw_stp.Tx_ID = 'sn:10447354119600060d001800cf281e583b';
hw_stp.Rx_ID = 'sn:104473dc599300131200210082672a4170';
hw_stp.N = 40e3;

%gnerate preamble
Preamble_1 = mseq(2,7);

Preamble_2 = mseq(2,7);
%generatiing source bits
source_bits = (randn(1,512*2))>0;

myModulator = comm.PSKModulator(...
        'ModulationOrder',  4,...  
        'PhaseOffset',      pi/4, ...     
        'BitInput',         true);         
 mod_symbols = myModulator(transpose(source_bits));

 %Alternatively
% data_symbolsQPSK = randi([0 mod_stp.M-1], 512, 1);
% mod_symbols = qammod(data_symbolsQPSK, mod_stp.M ,UnitAveragePower=true);

%mod_symbols = transpose(myModulator(transpose(source_bits)));
guard_bits = zeros(10,1);

tx_frame = [Preamble_1; Preamble_2; mod_symbols; guard_bits];

 gauss_filter = gaussdesign(0.3,flt_stp.span,flt_stp.sps);

tx_frame_up = upsample(tx_frame,flt_stp.sps);

filteredFrame_3 = conv(gauss_filter, tx_frame_up);


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
 Tx.transmitRepeat(filteredFrame_3); 
 %[frm,symbs_tx] = transmit(Tx,mod_stp,flt_stp);
 [smpls_rx_3,prmb_rx_3,smpls_dt_3,symbs_rx_3] = receive(Rx,mod_stp, flt_stp,hw_stp);


 %BT = 05

 gauss_filter = gaussdesign(0.5,flt_stp.span,flt_stp.sps);

tx_frame_up = upsample(tx_frame,flt_stp.sps);

filteredFrame_5 = conv(gauss_filter, tx_frame_up);


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
 Tx.transmitRepeat(filteredFrame_5); 
 %[frm,symbs_tx] = transmit(Tx,mod_stp,flt_stp);
 [smpls_rx_5,prmb_rx_5,smpls_dt_5,symbs_rx_5] = receive(Rx,mod_stp, flt_stp,hw_stp);

  %BT = 07
 gauss_filter = gaussdesign(0.7,flt_stp.span,flt_stp.sps);

tx_frame_up = upsample(tx_frame,flt_stp.sps);

filteredFrame_7 = conv(gauss_filter, tx_frame_up);


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
 Tx.transmitRepeat(filteredFrame_7); 
 %[frm,symbs_tx] = transmit(Tx,mod_stp,flt_stp);
 [smpls_rx_7,prmb_rx_7,smpls_dt_7,symbs_rx_7] = receive(Rx,mod_stp, flt_stp,hw_stp);

 %SNR USING EVM

 evm = comm.EVM;
 rmsEVM1 = evm(mod_symbols, symbs_rx/rms(symbs_rx));

 SNR = -20*log10(rmsEVM1/100);


%plotting power spectrum
figure()
[Pxx_tx3,F_tx3] = pwelch(smpls_rx_3,[],[],[],hw_stp.R,'centered');
plot(F_tx3,10*log10(Pxx_tx3));
grid on
hold on
[Pxx_tx5,F_tx5] = pwelch(smpls_rx_5,[],[],[],hw_stp.R,'centered');
plot(F_tx5,10*log10(Pxx_tx5));
hold on
[Pxx_tx7,F_tx7] = pwelch(smpls_rx_7,[],[],[],hw_stp.R,'centered');
plot(F_tx7,10*log10(Pxx_tx7));
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Rx signal')
legend('0.3','0.5','0,7')

%99 OBW
figure()
obw(smpls_rx,hw_stp.R,[],99);

%CCDF
figure()
ccdf = comm.CCDF(...
    'AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true);
[ccdfy,ccdfx,avg,peak] = ccdf([smpls_rx smpls_rx]);

pm = powermeter(ComputeCCDF=true);
averagePower = pm(smpls_rx); 
prob = probability(pm,3);
plotCCDF(pm,GaussianReference=true);
 title('CCDF Measurements QPSK')

plot(ccdf)
legend('QPSK')

%constellation diagram
scatterplot(mod_symbols);

scatterplot(symbs_rx/rms(symbs_rx));

%eye diagrams
eyediagram(smpls_rx, 64)

release(Tx);
release(Rx);