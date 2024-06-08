%Building a frame
% preamble
Pr = [mseq(2,7); mseq(2,7)];
Pr=Pr /rms(Pr);
%Symbols QPSK
sps = 8;
M = 16;
x = randi([0 M-1],512,1);
y = qammod(x,M);

qpsk_symbols = transpose(y);
qpsk_symbols=qpsk_symbols/rms(qpsk_symbols);
scatterplot(qpsk_symbols)
% the frame before filtering 
frm_qpsk = [zeros(128,1);Pr; qpsk_symbols';zeros(128,1)];
scatterplot(frm_qpsk)
% upsmapling

frm_qpsk_up = upsample(frm_qpsk,sps);

% NOw lets filter
alpha = 0.5;

span = 16;
fltr_rc = rcosdesign(alpha,span,sps,"sqrt");

frm= conv(fltr_rc,frm_qpsk_up);

[p, f]=pwelch(frm_samples_qpsk_rc,[],[],[],1e6,"centered");
plot(f,10*log(abs(p)));

%% Setting up Adalm Pluto
% 1. initialize two objects 
hw_stp.Tx_G=-10;
hw_stp.Rx_G=20;
hw_stp.F=2.4e9;
hw_stp.R=1e6;
% sd=findPlutoRadio
hw_stp.Tx_ID='sn:1044730a1997001610002a0036067c6324';
hw_stp.Rx_ID='sn:1044734c960500072400190048ae4f7dcb';
hw_stp.N=30e3; % they said 30e3;

mod_stp.T='QAM';
mod_stp.M=16;
mod_stp.N=512;

%The filter setup is defined differently for each part
flt_stp.T='RRC';
flt_stp.sps=8;
flt_stp.span=16;
flt_stp.alpha=0.5;
flt_stp.BT=0.3;
flt_stp.mtchd= true;


% Initialize Tx object with the given parameters
Tx = sdrtx( 'Pluto','Gain', hw_stp.Tx_G, ...
    'CenterFrequency', hw_stp.F, ...
    'BasebandSampleRate',hw_stp.R, ...
    'RadioID', hw_stp.Tx_ID);




%RxRadioID = PlutoRadioSerialNumbers(2);

% Initialize Rx object with the given parameters
Rx = sdrrx( 'Pluto', 'SamplesPerFrame', hw_stp.N, ...
    'OutputDataType', 'double', ...
    'CenterFrequency', hw_stp.F, ...
    'BasebandSampleRate',hw_stp.R,...
    'GainSource', 'AGC Fast Attack', ...
    'Gain', hw_stp.Rx_G, ...
    'RadioID', hw_stp.Rx_ID);



Tx.transmitRepeat(frm);
sig_rx = Rx();
release(Tx)
release(Rx)

figure()
plot(20*log10(abs(sig_rx))); % time domain
selected_frame = sig_rx(7901:16378);

[p, f]=pwelch(sig_rx,[],[],[],hw_stp.R,"centered"); %freqyebcy domain
plot(f,10*log10(p));

plot(20*log10(abs(selected_frame)))
 %energy detection
edge_crs=dtct_edge2(selected_frame)
 if 1
plot(20*log10(abs(selected_frame)));
hold on
plot(edge_crs,20*log10(abs(selected_frame)),'r*');
 end

 %% for match filtering
 fltr_rc = rcosdesign(alpha,span,sps,"sqrt");
 samples_mtchd = conv(fltr_rc,selected_frame);
 eyediagram(samples_mtchd,32)
 






% sig_rx_mod16 = sig_rx_mod4:
% sig_rx = sig_rx_mod4;
% plot(abs(sig_rx_mod16_withguard))
% 



one_frame = sig_rx(7431:15794);
plot(abs(one_frame))

k = length(one_frame);
der =[];
for i=2:k
 der(i) = one_frame(i)-one_frame(i-1);
end

plot(abs(der))
frm1 = conv(fltr_rc,upsample(Pr,sps));
g = xcorr(one_frame,frm1);
plot(abs(g));

dtct_edge2(one_frame)


plot(10*log10(abs(sig_rx)));

samples = one_frame(1280:end);
plot(abs(samples(1:sps*254)))

samples_mtchd = conv(fltr_rc,samples);
eyediagram(samples_mtchd(8*254:end),sps*4);

plot(abs(samples_mtchd))

symbols = downsample(samples_mtchd,sps,);

scatterplot(symbols)


prs = samples(51:1880);


Xf = fftshift(fft(prs.^2));
fs = hw_stp.R;
F = linspace(-fs/2,fs/2,length(Xf));
plot(F,abs(Xf))

foff = 13395.3/2;

% compensating for frequency offset
samples_mtchd_compenat = [];
for i=1:length(samples_mtchd)
    samples_mtchd_compenat(i)=samples_mtchd(i)*exp(j*2*pi*foff*i/fs);

end
eyediagram(samples_mtchd_compenat(sps*254:end),4*sps)