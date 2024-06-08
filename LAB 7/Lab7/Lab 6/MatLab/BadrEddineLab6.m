%% PATH-LOSS EXPONENT
% 2. build a baseband signal to transmit a single tone
sps = 8;
num_of_sym = 100;
sym = (ones(1,num_of_sym)+j*ones(1,num_of_sym))/sqrt(2);
temp      = repmat(sym',sps,1);
frm_samples_rect = reshape(temp,1,num_of_sym*sps);
frm_samples_rect = frm_samples_rect';

%3. Received signal power
Tx.transmitRepeat(frm_samples_rect); 
single_tone_rx1_los_d1 = Rx();
release(Tx)
release(Rx)

[p, f]=pwelch(single_tone_rx1_los_d1,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));
p1 = -155.966;
p2 = -161.224;
%4. non-LoS
p3 = -166.985;
p4 = -187.693;
%5. Comment on the results obtained in (3) and (4). How do you compare the path gains
%in the two cases?
%answer : the path gain with los is higher than that of the nlos
%the gain for los over nlos at d1 is -155.966-(-166.985)=11.019dB
%the gain for los over nlos at d2 is -161.224+187.693= 26.469dB
% the average is 18.7440
% Hence the gain with los is on average 18.744 dB greater than nlos


% 6. For a frequency of 915MHz
%6.3
[p, f]=pwelch(single_tone_rx2_nlos_d2,[],[],[],hw_stp.R,"centered");
plot(f,10*log(abs(p)));

p1 = -103.761;
p2 = -117.646;
%6.4. non-LoS
p3 = -124.016;
p4 = -132.419;
%6.5 comments
% we notice that with decreasing the frequency, the path loss is lower
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

%9. sendind and receiving at different sample rates
%9. Power Spectral Density
hw_stp.R = 15e6;
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
scatterplot(symbs_rx_rc14/h);
pr_tx = [Pr;Pr];
pr_rx = prmb_rx_rc14;

h = mean(pr_tx./pr_rx);

scatterplot(symbs_rx_rc14);
equalized_symbs = symbs_rx_rc14/h;
scatterplot(equalized_symbs);


recov_eq = zeros(512,1);
for ii=1:length(symbs_rx_rc14)

    if real(equalized_symbs(ii)<0)
        if(imag(equalized_symbs(ii))<0)
            recov_eq(ii) = exp(-j*3*pi/4);
        else
            recov_eq(ii) = exp(j*3*pi/4);
        end
    else
        if(imag(equalized_symbs(ii))<0)
            recov_eq(ii) = exp(-j*pi/4);
        else
            recov_eq(ii) = exp(j*pi/4);
    end
end
end

% Comparing the recovered symbols with the received symbols
g1 = xcorr(recov_eq,qpsk_symbols');
plot(abs(g1))

cm = recov_eq-qpsk_symbols';
max(abs(cm))
scatterplot(cm);
scatterplot(recov_eq);

