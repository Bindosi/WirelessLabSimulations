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

%qpsk modulation
mod_symbols = pskmod(data, mod_stp.M );

% Upsampling
%tx_frame_up = upsample(mod_symbols,flt_stp.sps);

%creatıng our fılter
filter = rcosdesign(flt_stp.alpha, flt_stp.span, flt_stp.sps,"normal");

%pefroming convolution to get the pulse shaped frame
filteredFrame = conv(filter, mod_symbols);

%creating time interval for the phase shifts
time_interv = 0:(1/flt_stp.R): (((length(mod_symbols)-1)/flt_stp.R));

%downsampling to get original number of samples with offset sample of 1
tx_frame_down = filteredFrame(1+flt_stp.span*flt_stp.sps/2-1: ...
    length(filteredFrame)-(flt_stp.sps*flt_stp.span/2)-1);

%introducing frequency offset
offset_tx_frame_down = tx_frame_down.*exp(1i*2*pi*foffset.*time_interv');

%fft shift to get frequency components
sampx = fftshift(fft(offset_tx_frame_down.^mod_stp.M));

scatterplot(mod_symbols);
scatterplot(offset_tx_frame_down)

delta_f = (flt_stp.R)/(mod_stp.N);

freq_interv=(-(flt_stp.R)/2):delta_f:((flt_stp.R)/2)-1;

figure()
grid on
plot(freq_interv,abs(sampx));
s_froff  =    (sampx == max(sampx))/mod_stp.M;

