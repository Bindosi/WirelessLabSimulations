%% PATH-LOSS EXPONENT

% 1. initialize two objects 
hw_stp.Tx_G=-10;
hw_stp.Rx_G=20;
hw_stp.F=2.4e9;
hw_stp.R=1e6;
% sd=findPlutoRadio
hw_stp.Tx_ID='sn:10447354119600060b003a0007241ed971';
hw_stp.Rx_ID='sn:1044735411960004f9ff0c00f1ff1efda3';
hw_stp.N=30e3; % they said 30e3;

mod_stp.T='QPSK';
mod_stp.M=4;
mod_stp.N=512;

%The filter setup is defined differently for each part
flt_stp.T='RRC';
flt_stp.sps=8;
flt_stp.span=12;
flt_stp.alpha=0.1;
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

% 2. build a baseband signal to transmit a single tone
sps = 8;
num_of_sym = 100;
sym = (ones(1,num_of_sym)+j*ones(1,num_of_sym))/sqrt(2);
temp      = repmat(sym',sps,1);
frm_samples_rect = reshape(temp,1,num_of_sym*sps);
frm_samples_rect = frm_samples_rect';

Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx = Rx();
release(Tx)
release(Rx)
figure()
[p, f]=pwelch(single_tone_rx,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

%3. Received signal power
Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx1_los_d1 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx1_los_d1,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx1_los_d2 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx1_los_d2,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));
%4. non-LoS
Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx1_nlos_d1 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx1_nlos_d1,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx1_nlos_d2 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx1_nlos_d1,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

%6. Repeating 3 and 4 for a different frequency
% change center frequency then transmit and capture

Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx2_los_d1 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx2_los_d1,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx2_los_d2 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx2_los_d2,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));
%4. non-LoS
Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx2_nlos_d1 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx2_nlos_d1,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx2_nlos_d2 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx2_nlos_d2,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

%% FREQUENCY-SELECTIVITY
% set the center frequncy back, use ‘AGC Fast Attack’

% 8. Building a frame
% preamble
Pr = mseq(2,7);
%Symbols QPSK
source_bits = (randn(1,512))>0;
bpsk_symbols = source_bits*2-1;

source_bits = (randn(1,512))>0;
myModulator = comm.PSKModulator(...
        'ModulationOrder',  2,...  
        'PhaseOffset',      0, ...     
        'BitInput',         true);

bpsk_symbols = transpose(myModulator(transpose(source_bits)));

scatterplot(bpsk_symbols);


% the frame before filtering 
frm_bpsk = [Pr; Pr; bpsk_symbols'];

% upsmapling

frm_bpsk_up = upsample(frm_bpsk,sps);

% NOw lets filter
alpha = 0.1;
span = 12;
sps = 8;
fltr_rc = rcosdesign(alpha,span,sps,"sqrt");
frm_samples_bpsk_rc = conv(fltr_rc,frm_bpsk_up);

% Adding Guard samples
frm_samples_bpsk_rc = [zeros(200,1); frm_samples_bpsk_rc ; zeros(200,1)];
%2. send the frame with different sample-rates (e.g., 0.5Msps, 1Msps, 10Msps, 15Msps).

% first sample rate
frm_samples_bpsk_rc = frm_samples_bpsk_rc;

Tx.transmitRepeat(frm_samples_bpsk_rc); 

[smpls_rx_rc91,prmb_rx_rc91,smpls_dt_rc91,symbs_rx_rc91]=receive(Rx,mod_stp,flt_stp,hw_stp);
release(Tx)
release(Rx)

scatterplot(symbs_rx_rc91);
[p, f]=pwelch(smpls_rx_rc91,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

% second sample rate
Tx.transmitRepeat(frm_samples_bpsk_rc); 
[smpls_rx_rc92,prmb_rx_rc92,smpls_dt_rc92,symbs_rx_rc92]=receive(Rx,mod_stp,flt_stp,hw_stp);
release(Tx)
release(Rx)

[p, f]=pwelch(smpls_rx_rc92,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

% third sample rate
Tx.transmitRepeat(frm_samples_bpsk_rc); 
[smpls_rx_rc93,prmb_rx_rc93,smpls_dt_rc93,symbs_rx_rc93]=receive(Rx,mod_stp,flt_stp,hw_stp);
release(Tx)
release(Rx)

[p, f]=pwelch(smpls_rx_rc93,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

%fourth sample rate
Tx.transmitRepeat(frm_samples_bpsk_rc); 
[smpls_rx_rc94,prmb_rx_rc94,smpls_dt_rc94,symbs_rx_rc94]=receive(Rx,mod_stp,flt_stp,hw_stp);
release(Tx)
release(Rx)

[p, f]=pwelch(smpls_rx_rc94,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

%% Channel Estimation
% sending and capturing qpsk data

%Building a frame
% preamble
Pr = mseq(2,7);
%Symbols QPSK
source_bits = (randn(1,2*512))>0;
myModulator = comm.PSKModulator(...
        'ModulationOrder',  4,...  
        'PhaseOffset',      pi/4, ...     
        'BitInput',         true);

qpsk_symbols = transpose(myModulator(transpose(source_bits)));


% the frame before filtering 
frm_qpsk = [Pr; Pr; qpsk_symbols'];

% upsmapling

frm_qpsk_up = upsample(frm_qpsk,sps);

% NOw lets filter
alpha = 0.1;
span = 12;
fltr_rc = rcosdesign(alpha,span,sps,"sqrt");

frm_samples_qpsk_rc = conv(fltr_rc,frm_qpsk_up);
frm_samples_qpsk_rc = [zeros(200,1); frm_samples_qpsk_rc ; zeros(200,1)];

Tx.transmitRepeat(frm_samples_qpsk_rc); 
[smpls_rx_rc14,prmb_rx_rc14,smpls_dt_rc14,symbs_rx_rc14]=receive(Rx,mod_stp,flt_stp,hw_stp);
release(Tx)
release(Rx)
scatterplot(symbs_rx_rc14);
pr_tx = [Pr;Pr];
pr_rx = prmb_rx_rc14;

h = [];
h(1) = pr_rx(1)/pr_tx(1);
M = 10;
for n==2:M
    s = 0;
    for i==1:n
        s = s+pr_tx(i)*h(n-i);
    end
    h(n)= (pr_rx(n)-s)/pr_tx(1);
end

eq(0) = symbs_rx_rc14(1);

%% Questions


%% PATH-LOSS EXPONENT
