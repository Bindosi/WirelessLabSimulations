clear; clc;

% Edge detection parameters
sys.fs = 1e6;
sys.N = 512;
sys.M = 4;
sys.M_prmbl = 2;
sys.sps = 8;
sys.alpha = 0.5;
sys.span = 16;
load('OFDM_RVC_CPSize4.mat')

m_sq        = mseq(2,7);
prmbl       = [m_sq; m_sq];
prmbl       = prmbl/rms(prmbl);




% time domain
figure()
plot((abs(Y)));
figure()

%frequency domain plot of the received singals
[Pxx_tx_q_a,F_tx_q_a] = pwelch(Y,[],[],[],sys.fs,'centered');
plot(F_tx_q_a,10*log10(Pxx_tx_q_a));
title('Power Spectrum of the Received Signal')
grid on

% Selecting the frame
selected_frame = Y(3132:5021);

%% Edge detection and Testing the Edge Detection Process
coarse_edge_det     = edge_detect(selected_frame,sys);
if (1) % Edge detection result
        figure()
        plot(20*log10(abs(selected_frame)));
        hold on
        plot(coarse_edge_det , 20*log10(abs(selected_frame)),'r*');
        hold off
end

 %% Estimation of frequency offset and compensation
    symbls_off = 10;

    indx_start = coarse_edge_det + symbls_off * sys.sps;
    indx_end = indx_start + (length(prmbl)-symbls_off)*sys.sps;

    preambles = selected_frame(indx_start:indx_end);

    coarse_offset = est_coarse_fre(preambles, sys);
    % Compensation of frequency offset
    t = (1:length(selected_frame))'/sys.fs;
    % compensation of the frequency offset
    smpls_rx_q = selected_frame.*exp(-2i*pi*t*coarse_offset);
         % scatterplot(selected_frame)
         % scatterplot(smpls_rx_q)
         % eyediagram(smpls_rx_q,128)
%% Match filtering
 fltr = rcosdesign(sys.alpha, sys.span, sys.sps, 'sqrt');
 sampls_rx_q_mt = conv(selected_frame, fltr);
 % eyediagram(sampls_rx_q_mt,256)
 % scatterplot(sampls_rx_q_mt)
%% Time Synchronization (Downsampling and fine edge detection)

% (0) Find the frame edges with guard symbols
% if(1)
% symbols_off             = 20;
% frame_size              = length(prmbl)+sys.N;
% index_start             = coarse_edge_det - symbols_off*sys.sps;
% index_end               = coarse_edge_det + frame_size*sys.sps + symbols_off*sys.sps;
% singal_downsampled      = sampls_rx_q_mt(index_start:index_end);
% 
% % (2) Downsample
% [prmbl_rx, symbols_rx]  = dwn_smpl2(singal_downsampled, prmbl,symbols_off,sys);
% frame_rx                =   [prmbl_rx; symbols_rx];
%       scatterplot(prmbl_rx)
%       scatterplot(symbols_rx)
% end

% %% Downsampling without frequecy offset 
% (0) Find the frame edges with guard symbols
if(1)
symbols_off             = 20;
frame_size              = length(prmbl)+sys.N;
index_start             = coarse_edge_det - symbols_off*sys.sps;
index_end               = coarse_edge_det + frame_size*sys.sps + symbols_off*sys.sps;
singal_downsampled      = sampls_rx_q_mt(index_start:index_end);

% (2) Downsample
[prmbl_rx, symbols_rx]  = dwn_smpl2(singal_downsampled, prmbl,symbols_off,sys);
frame_rx                =   [prmbl_rx; symbols_rx];
%  off_ifft  = est_frq2(prmbl_rx, prmbl,sys);
% 
% ts = (1:length(symbols_rx))'/sys.fs;
% tp = (1:length(prmbl_rx))'/sys.fs;
% 
%  symbols_rx = symbols_rx.*exp(-2i*pi*ts*off_ifft);
%  prmbl_rx   = prmbl_rx.*exp(-2i*pi*tp*off_ifft);
     scatterplot(sampls_rx_q_mt) 
     scatterplot(prmbl_rx)
     scatterplot(symbols_rx)
end
 %% Phase Ambiguity/Channel Estimation
 estimated_h         = mean(prmbl_rx./prmbl);    
    % h_cap               = prmbl_rx.* conj(prmbl)/abs(prmbl).^2;
    % estimatedChannel    = conv(prmb_rx_q14, fliplr(Preamble_q14)); % Convolution Approach
    %channel_check       = isequal(estimated_h,channel);

% Equalization of received symbols to get the original symbols
eq_smpls_rx         =   symbols_rx.* conj(estimated_h(1));
estimated_pr        =   prmbl_rx.* conj(estimated_h(1));

phase_ambiguity     =   angle(estimated_h)/pi;
phase_am_magnitude =    abs(estimated_h)
 scatterplot(eq_smpls_rx) 
 scatterplot(estimated_pr)


% estimating bit error rate
for k = 1:length(prmbl)
    BER(k) = (length(find(prmbl ~=estimated_pr)))/sys.fs;
end

% Symbol detection (nearest neighbor)
M = 16;
% Demodulate the received signal
demodulatedData = qamdemod(eq_smpls_rx, M, 'UnitAveragePower', true);
detectedSymbols = zeros(size(demodulatedData));
for i = 1:length(demodulatedData)
    distances = abs(demodulatedData(i) - qammod(0:M-1, M, 'UnitAveragePower', true));
    [~, index] = min(distances);
    detectedSymbols(i) = index - 1; % Adjust to 0-based indexing
end
%claculating EVM and SNR
EVM = sqrt(mean(abs(eq_smpls_rx-detectedSymbols).^2));

SNR = -20*log10(EVM/100);

%% Edge detectıon
function edge_index  = edge_detect(selected_frame,sys)
    window_size         = sys.sps*16;               %detection window
    absolute_values     = abs(selected_frame).^2;
    vector2matrix       = vec2mat(absolute_values, window_size);
    derivatives         = diff(sum(vector2matrix.'));
    plot(derivatives);
     [~, indx] = max(derivatives);
     plot(derivatives)
     edge_index = indx*window_size;
end

%% Frequency estimation 
function off = est_coarse_fre(preambles,sys)
    
    sampx = fftshift(fft(preambles.^2,40*length(preambles)));
    delta_f = (sys.fs)/length(sampx);

    freq_interv=(-(sys.fs))/2:delta_f:((sys.fs)/2)-1;
    figure()
    plot(freq_interv/2, abs(sampx));
    off  = freq_interv(sampx==max(sampx))/2;
end

%% FFT based frequency offset estimation done at symbol rate, no compensation
function off_fft = est_frq2(rec_preambles, prmbl, sys)
fk = sys.fs/sys.sps; 
    % step 1: coarse estimate
    symbs_off_pwr = rec_preambles.*prmbl;
    signal_fft = fftshift(fft(symbs_off_pwr));
    kk = find(abs(signal_fft) == max(abs(signal_fft)));
    Nk = length(signal_fft);
    
    % samples
    y1 = (signal_fft(kk+1));
    y0 = (signal_fft(kk+0));
    ym1 = (signal_fft(kk-1));
    % ratios
    a1 = real(y1/y0);
    am1 = real(ym1/y0);
    
    %indicators
    sg1 = -a1/(1-a1);
    sgn1 = am1/(1-am1);
    
    if sg1 > 0 && sgn1>0
            sg = sg1;
        else
            sg = sgn1;
    end
    df = fk/Nk;
    freq = -fk/2:df:fk/2;
    freq = freq(1:end-1);
    off_fft = freq(kk);             % 
    off_fft = (off_fft+df*sg);      % Fractional Frequency Offset
end
%% downampling function
function [prmbl_rx, symbls_rx] = dwn_smpl2(smpls, prmbl, symbols_off, sys) %%Power-based
    sps = sys.sps;
    indices = zeros(sps,1);
    for ii = 1:sps
        crr = sum(abs(smpls(ii:sps:end)).^2);
        indices(ii) = crr;
    end
    opt = find(indices ==  max(indices));
    symbs = downsample(smpls(opt:end),sps);

    indx = 1;
    prmbl_length = length(prmbl);
    crr_prv = 0;
    for i = 1:2*symbols_off
        crr = abs(prmbl'*symbs(i:i+prmbl_length-1));
        correl(i) = crr; 
        if crr > crr_prv
            crr_prv = crr;
            indx = i;
        end
        
    end
    figure()
    plot(correl);
    title('Correlation for the preambles')
    indx_start  = indx+prmbl_length;
    symbls_rx   = symbs(indx_start:indx_start+sys.N-1);
    prmbl_rx    = symbs(indx:indx+prmbl_length-1);
end

