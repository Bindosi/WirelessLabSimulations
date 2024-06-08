mod_stp.T = 'QPSK';
mod_stp.M = 4;
mod_stp.N = 512;


flt_stp.T = 'RC';
flt_stp.sps = 8;
flt_stp.span = 12;
flt_stp.alpha = 0.9;

%generate symbols
data =  randi([0 mod_stp.M-1],mod_stp.N,1);

%qpsk modulation
mod_symbols = pskmod(data, mod_stp.M, pi/mod_stp.M );

% Upsampling
tx_frame_up = upsample(mod_symbols,flt_stp.sps);

%creatıng our fılter
filter = rcosdesign(flt_stp.alpha, flt_stp.span, flt_stp.sps,"normal");

%pefroming convolution to get the filtered frame
filteredFrame = conv(filter, tx_frame_up);

% tx_frame_down = downsample(filteredFrame(1+ (flt_stp.span*flt_stp.sps/2):...
    % flt_stp.sps: length(filteredFrame) - 1),flt_stp.sps,1);

tx_frame_down = downsample(filteredFrame(flt_stp.span*flt_stp.sps/2-1: ...
    length(filteredFrame)-(flt_stp.sps*flt_stp.span/2+1)) ,flt_stp.sps,0);

%Correlation
[c,lags] = xcorr(tx_frame_down, mod_symbols);
plot(lags,c);

% constelation
scatterplot(mod_symbols);
scatterplot(tx_frame_down);

