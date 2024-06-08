%Mod 4QAM
% Filter Setup
s.T ="RC";
s.sps = 8;
s.span = 12;
s.alpha = 0.5;
%Preamble
Pr1 = mseq(2,7);
Pr2 = mseq(2,7);
%QAM
m = 4;
data_symbols = randi([0 m-1], 512, 1);
mod_symbols = qammod(data_symbols, m);

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
temp4qm = [Pr1; Pr2; mod_symbols];
frm4qmTx = fltr(temp4qm, s);
mod_stp.T='QAM';
mod_stp.M=4;
mod_stp.N=512;
Tx.transmitRepeat(frm4qmTx); 
[smpls_rx4qm,prmb_rx4qm,smpls_dt4qm,symbs_rx4qm]=receive(Rx,mod_stp,s,hw_stp);

scatterplot(symbs_rx4qm);

% 16-QAM
% Filter Setup
s.T ="RC";
s.sps = 8;
s.span = 12;
s.alpha = 0.5;
%Preamble
Pr1 = mseq(2,7);
Pr2 = mseq(2,7);
%QAM
m = 16;
data_symbols16 = randi([0 m-1], 512, 1);
mod_symbols16 = qammod(data_symbols16, m ,UnitAveragePower=true);

scatterplot(mod_symbols16)

hw_stp.Tx_G=0;
hw_stp.Rx_G=40;
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
    'GainSource', 'AGC Fast Attack', ...
    'RadioID', hw_stp.Rx_ID);
temp16qm = [Pr1; Pr2; mod_symbols16];
frm16qmTx = fltr(temp16qm, s);
mod_stp.T='QAM';
mod_stp.M=m;
mod_stp.N=512;
Tx.transmitRepeat(frm16qmTx); 
[smpls_rx16qm,prmb_rx16qm,smpls_dt16qm,symbs_rx16qm]=receive(Rx,mod_stp,s,hw_stp);
scatterplot(symbs_rx16qm);

% 64-QAM
% Filter Setup
s.T ="RC";
s.sps = 8;
s.span = 12;
s.alpha = 0.5;
%Preamble
Pr1 = mseq(2,7);
Pr2 = mseq(2,7);
%QAM
m = 64;
data_symbols64 = randi([0 m-1], 512, 1);
mod_symbols64 = qammod(data_symbols64, m,UnitAveragePower=true);

hw_stp.Tx_G=0;
hw_stp.Rx_G=30;
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
temp64qm = [Pr1; Pr2; mod_symbols64];
frm64qmTx = fltr(temp64qm, s);
figure
plot(abs(frm64qmTx))
mod_stp.T='QAM';
mod_stp.M=m;
mod_stp.N=512;
Tx.transmitRepeat(frm64qmTx); 
[smpls_rx64qm,prmb_rx64qm,smpls_dt64qm,symbs_rx64qm]=receive(Rx,mod_stp,s,hw_stp);
figure
plot(abs(smpls_rx64qm))
scatterplot(symbs_rx64qm);
%% Answering Questions.
% Q2 Power spetral Densities

figure()
[Pxx_tx4qm,F_tx4qm] = pwelch(smpls_rx4qm,[],[],[],hw_stp.R,'centered');
plot(F_tx4qm,10*log10(Pxx_tx4qm));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)');
title('Power Spectrum of the 4QAM signal');
obw(smpls_rx4qm,hw_stp.R,[],99);
[P4qm, F4qm] = pspectrum(smpls_rx4qm, hw_stp.R);
plot(F4qm,10*log10(P4qm));
figure()
[Pxx_tx16qm,F_tx16qm] = pwelch(smpls_rx16qm,[],[],[],hw_stp.R,'centered');
plot(F_tx16qm,10*log10(Pxx_tx16qm));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)');
title('Power Spectrum of the 16QAM signal');
figure()
obw(smpls_rx16qm,hw_stp.R,[],99);

figure()
obw(smpls_rx4qm,hw_stp.R,[],99);

figure()
 obw(smpls_rx64qm,hw_stp.R,[],99);


figure()
[Pxx_tx64qm,F_tx64qm] = pwelch(smpls_rx64qm,[],[],[],hw_stp.R,'centered');
plot(F_tx64qm,10*log10(Pxx_tx64qm));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)');
title('Power Spectrum of the 64QAM signal');


% USing pspectrum
figure()
[P4qm, F4qm] = pspectrum(smpls_rx4qm, hw_stp.R);
plot(F4qm,10*log10(P4qm));
title('Pspectrum of the 4QAM signal');

figure()
[P16qm, F16qm] = pspectrum(smpls_rx16qm, hw_stp.R);
plot(F16qm,10*log10(P16qm));
title('Pspectrum of the 16QAM signal');

figure()
[P64qm, F64qm] = pspectrum(smpls_rx64qm, hw_stp.R);
plot(F64qm,10*log10(P64qm));
title('Pspectrum of the 64QAM signal');

e = comm.EVM
e = comm.EVM
evm4 = e(temp4qm(end-511:end),symbs_rx4qm/rms(symbs_rx4qm));

e = comm.EVM
e = comm.EVM
evm16 = e(temp16qm(end-511:end),symbs_rx16qm/rms(symbs_rx16qm));

e = comm.EVM
e = comm.EVM
evm64 = e(temp64qm(end-511:end),symbs_rx64qm/rms(symbs_rx64qm));
% Q3 contellation and Eye Diagrams

scatterplot(symbs_rx4qm);
title('Constellation diagram of the 4QAM signal');
scatterplot(symbs_rx16qm);
title('Constellation diagram  of the 16QAM signal');
scatterplot(symbs_rx64qm);
title('Constellation diagram  of the 64QAM signal');

eyediagram(symbs_rx4qm, 4);
title('Eye diagram of the 4QAM signal');
eyediagram(symbs_rx16qm, 4);
title('Eye diagram of the 16QAM signal');
eyediagram(symbs_rx64qm, 4);
title('Eye diagram of the 64QAM signal');

v=200;
eyediagram(smpls_rx16qm(end-300-v:end-v), 8*4);
eyediagram(smpls_rx64qm, 8*10);
% Becaues the receiver undoes the filtering, it actually gives better
% samples than the transmitted samples after the filter.
eyediagram(frm4qmTx(end-512:end), 8*4);