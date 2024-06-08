
clc, clear;
% code for lab 7
mod_stp.T = '16QAM';
mod_stp.M = 16;
mod_stp.N = 512;

% Edge detection parameters
sys.fs = 1e6;
sys.N = 512;
sys.M = 16;
sys.M_prmbl = 2;
sys.sps = 8;
sys.alpha = 0.5;
sys.span = 16;

% Filter set up
flt_stp.T = 'RRC';
flt_stp.sps = 8;
flt_stp.span = 16;
flt_stp.alpha = 0.5;
flt_stp.BT = 0.7;
flt_stp.mtchd = true;

%Hardware Setup
hw_stp.Tx_G  = -10;
hw_stp.Rx_G  = 20;
hw_stp.F     = 2.35e9;
hw_stp.R     = 0.5e6;
hw_stp.Tx_ID ='sn:10447354119600060d001800cf281e583b';
hw_stp.Rx_ID ='sn:104473dc599300131200210082672a4170';
hw_stp.N = 30e3;
%% Setting up Transmitter and Receiver
% setting up transmit object
        Tx  = sdrtx('Pluto', 'Gain',hw_stp.Tx_G,...
            'CenterFrequency',hw_stp.F,...
            'BasebandSampleRate',hw_stp.R,...
            'RadioID',hw_stp.Tx_ID);

% setting up receive object
        Rx = sdrrx('Pluto', 'SamplesPerFrame', hw_stp.N,...
            'OutputDataType','double',...
            'CenterFrequency',hw_stp.F,...
            'BasebandSampleRate',hw_stp.R,...
            'GainSource','Manual',... %AGC Fast Attack
            'Gain', hw_stp.Rx_G,...
            'RadioID',hw_stp.Rx_ID);

%% Generating data and creating a frame
data        = randi([0 mod_stp.M-1], mod_stp.N,1);
symbs       = qammod(data, mod_stp.M);
symbs       = symbs/rms(symbs);
scatterplot(symbs);

% Preamble Generation
m_sq        = mseq(2,8);
prmbl       = [m_sq; m_sq];
prmbl       = prmbl/rms(prmbl);

% Forming the Frame
grd         = zeros(128,1);
frame       = [grd; prmbl; symbs; grd];
frame_up    = upsample(frame, flt_stp.sps);
fltr        = rcosdesign(flt_stp.alpha, flt_stp.span,flt_stp.sps, 'sqrt');
filtered_frame        = conv(fltr, frame_up);
    % t = 1:length(filtered_frame);
    % figure();
    % plot(filtered_frame(1600:1800));
%% transmitting and receiving 

%%%%%%%%% USING SIMULATION FOR TESTING THE CODE %%%%%%%%

channel             = rand(5,1)+1j*rand(5,1);
smpls_rx_q_a        = conv(channel, filtered_frame);
mutli_frames = [smpls_rx_q_a; smpls_rx_q_a; smpls_rx_q_a];

smpls_rx_q  = mutli_frames;

%%%%%%%%% USING SDR %%%%%%%

% % Sending to SDR
% Tx.transmitRepeat(filtered_frame); 
% 
% % Receiving for SDR
% [smpls_rx_q_a, prmb_rx_q_a, smpls_dt_q_a, symbs_rx_q_a] = receive(Rx,mod_stp, flt_stp,hw_stp);

%% observe the received signal in time and frequency domain

% Power spectrum in frequency domain
figure()
[Pxx_tx_q_a,F_tx_q_a] = pwelch(smpls_rx_q_a,[],[],[],hw_stp.R,'centered');
plot(F_tx_q_a,10*log10(Pxx_tx_q_a));
title('Power Spectrum of the Received Signal')
grid on

% Time domain of the received signal
t = 1: length(smpls_rx_q);

figure()
plot(20*log10(abs(smpls_rx_q)));

%% Coarse (Energy Based) Edge Detection
coarse_edge_det     = edge_detect(smpls_rx_q,sys);

%% Testing the Edge Detection Process
if (1)                      % Edge detection result
    plot(20*log10(abs(smpls_rx_q)));
    hold on
     plot(coarse_edge_det , 20*log10(abs(smpls_rx_q)),'r*');
     hold off
end
%% Estimation of frequency offset, normally done for  low end equipment
if (1) % Estimation of frequency offset
    symbls_off = 5;
    indx_start = coarse_edge_det+symbls_off+sys.sps;
    indx_end = indx_start + (length(prmbl)-symbls_off)*sys.sps;
    sgnal_off = smpls_rx_q(indx_start:indx_end);
    coarse_offset = est_coarse_fre(sgnal_off, sys);
    % Compensation of frequency offset
    t = (1:length(smpls_rx_q))'/sys.fs;
    smpls_rx_q = smpls_rx_q.*exp(-2i*pi*t*coarse_offset);
end
scatterplot(smpls_rx_q)
%% Match filtering, for high end equipment that is enough
fltr = rcosdesign(sys.alpha, sys.span, sys.sps, 'sqrt');
smpls_rx_q = conv(smpls_rx_q,fltr);

%% Freqency Offset Compensation (Duplicated m-sequence)
if 1
    off_duplicated = est_frq1(smpls_rx_q,coarse_edge_det,sys);
    t = (1:length(smpls_rx_q))'./sys.fs;
    smpls_rx_q = smpls_rx_q.*exp(-2i*pi*t*coarse_offset);
end

%% Time Synchronization

% (1) Find the frame edges with guard symbols
symbols_off             = 20;
frame_size              = length(prmbl)*sys.N;
index_start             = coarse_edge_det - symbols_off*sys.sps;
index_end               = coarse_edge_det + frame_size*sys.sps + symbols_off*sys.sps;
singal_downsampled      = smpls_rx_q(index_start:index_end);

% (2) Downsample
[pmble_rx, symbols_rx]  = dwn_smpl2(singal_downsampled,prmbl,symbols_off,sys);
frame_rx                =   [pmble_rx; symbols_rx];
scatterplot(prmbl_rx)
scatterplot(symbols_rx)
%% Frequency offset compensation FFT based (done at symbol rate)
    if 1
        off_ifft  = est_frq2(prmbl_rx, prmbl,sys);
         t = (1:length(frame_rx_q))'./sys.fs*sys.sps;
         frame_rx = frame_rx.*exp(-2i*pi*t*off_ifft);

         prmbl_rx = frame_rx(1:length(prmbl));
         symbls_rx = frame_rx(length(prmbl)+1:end);
    end

%% Phase ambiguity 
ang = fx_rot(prmbl_rx, prmbl);
prmbl_rx = prmbl_rx*exp(-1i*ang);
symbls_rx = symbls_rx*exp(-1i*ang);
scatterplot(symbls_rx);
scatterplot(prmbl_rx);
%% Energy Based Edge Detection Function
function edge_index  = edge_detect(smpls_rx_q,sys)
    window_size         = sys.sps*16;               %detection window
    absolute_values     = abs(smpls_rx_q).^2;
    vector2matrix       = vec2mat(absolute_values, window_size);
    derivatives         = diff(sum(vector2matrix.'));
    plot(derivatives);
     [~, indx] = max(derivatives);
     plot(derivatives)
     edge_index = indx*window_size;
end
%% Coarse frequency offset detection
function off = est_coarse_fre(smpls,sys)
    M = sys.M_prmbl;
    fs = sys.fs;
    sampls_pwr  = smpls.^M;  % this will make the signal a DC signal,                      
    sampls_fft  = abs(fftshift(fft(sampls_pwr,5*sampls_pwr )));  % FFT will give a peak at 0 frequency
    df          = fs/length(sampls_fft);
    freq        =   -fs/2:df:fs/2-1;
    plot(freq,sampls_fft);
    off         = freq(samples_fft == max(samples_fft))/M;

end

%% FFT based frequency offset estimation done at symbol rate, no compensation
function off_fft = est_frq2(symbs_off, prmbl, sys)
    fk = sys.fs/sys.sps;
    
    % step 1: coarse estimate
    symbs_off_pwr =symbs_off.*prmbl;
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
    sg1 = -a/(1-a1);
    sgn1 = am1/(1-am1);
    
    if sg1 > 0 && sgb1>0
            sg = sg1;
        else
            sg = sgn1;
    end
    df = fk/Nk;
    freq = -fk/2:df:fk/2;
    freq = freq(end-1);
    off_fft = freq(kk);             % 
    off_fft = (off_fft+df*sg);      % Fractional Frequency Offset
end


%% downampling function
function [prmbl_rx, symbls_rx] = dwn_smpl2(smpls, prmbl, off, sys) %%Power-based
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
    for i = 1:2*off
        crr = abs(prmbl'*symbs(i:i+prmbl_length-1));
        if crr > crr
            crr_prv = crr;
            indx = i;
        end
    end
    indx_start = indx+prmbl_length;
    symbls_rx = symbs(indx_start:indx_start+sys.N-1);
    prmbl_rx = symbs(indx:indx+prmbl_length-1);
end

function [prmbl_rx, symbls_rx] = dwn_smpl(smpls, prmbl, off, sys)
end
function off = est_frq1(sgnl_rx, prmbl, edg_crs, sys)
end
function ang = fx_rot(prmbl_rx, prmbl)
end








