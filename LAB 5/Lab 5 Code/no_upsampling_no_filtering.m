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

modData = real(pskmod(data,mod_stp.M ));

%creating time interval for the phase shifts
time_interv = 0:(1/flt_stp.R): (((length(mod_symbols)-1)/flt_stp.R));

%introducing frequency offset
offset_tx_frame = mod_symbols.*exp(1i*2*pi*foffset.*time_interv');

%fft shift to get frequency components
sampx = fftshift(fft(offset_tx_frame.^mod_stp.M));

%generating scatter plots
scatterplot(mod_symbols);           %for the original BPSK symbols
scatterplot(offset_tx_frame)        %for the received BPSK symbols

%determining delta f frequency spacing
delta_f = (flt_stp.R)/(mod_stp.N);  

%frequency interval for plotting diagram to see frequency offset
freq_interv=(-(flt_stp.R)/2):delta_f:((flt_stp.R)/2)-1;

figure()
grid on
plot(freq_interv,abs(sampx));
s_froff  =    (sampx == max(sampx))/mod_stp.M;

