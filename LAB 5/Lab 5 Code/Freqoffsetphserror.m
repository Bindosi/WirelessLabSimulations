clc 
clear all

mod_stp.T = 'BPSK';
mod_stp.M = 2;
mod_stp.N = 512;


flt_stp.T = 'RC';
flt_stp.sps = 8;
flt_stp.span = 12;
flt_stp.alpha = 0.9;
flt_stp.R     = 25e3;
foffset = 320;

%generate symbols
data =  randi([0 mod_stp.M-1],mod_stp.N,1);

%BPSK modulation
mod_symbols = pskmod(data, mod_stp.M );

% Alternatively you could get BPSK modulated symbols using
mod_symbols_alternate =  2*randi([0 mod_stp.M-1],mod_stp.N,1)-1;

% Upsampling
tx_frame_up = upsample(mod_symbols,flt_stp.sps);

%creatıng our fılter
filter = rcosdesign(flt_stp.alpha, flt_stp.span, flt_stp.sps,"normal");

%pefroming convolution to get the filtered frame
filteredFrame = conv(filter, tx_frame_up);

%creating time interval for the phase shift
time_interv = 0:(1/flt_stp.R): ((length(filteredFrame)-1)/flt_stp.R);


%introducing frequency offset
mod_phase_offset = filteredFrame.*exp(-1i*2*pi*foffset*time_interv)';

%downsampling the received signal
offset_tx_frame_down = downsample(mod_phase_offset(1+flt_stp.span*flt_stp.sps/2: ...
    length(mod_phase_offset)-(flt_stp.sps*flt_stp.span/2+1)) ,flt_stp.sps,0);

%fft shift to get frequency components
sampx = fftshift(fft(offset_tx_frame_down.^mod_stp.M));

scatterplot(mod_symbols);
scatterplot(offset_tx_frame_down)

delta_f = (flt_stp.R)/(mod_stp.N);

freq_interv=(-(flt_stp.R)/2):delta_f:((flt_stp.R)/2)-1;

figure()
plot(freq_interv/mod_stp.M, abs(sampx));
s_froff  = freq_interv(sampx==max(sampx))/mod_stp.M;

coarse_frequency_offset = abs(s_froff);