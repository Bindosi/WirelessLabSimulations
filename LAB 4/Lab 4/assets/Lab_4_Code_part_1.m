
% Filter Setup
s.T ="RC";
s.sps = 8;
s.span = 12;
s.alpha = 0.9;


%gnerate preamble
Preamble_1 = mseq(2,7);

Preamble_2 = mseq(2,7);
%generatiing source bits
source_bits = (randn(1,512))>0;

%qpsk modulation
    myModulator = comm.PSKModulator(...
        'ModulationOrder',  mod_stp.M,...   % here the second input determines how to
        'PhaseOffset',      pi/4, ...      % rotate (phase offset) the constellation
        'BitInput',         true);         % point in the complex plane. If we don't
    % rotate, it provides QPSK costellations
    % on the I and Q axis like [1 j -1 -j].
    mod_symbols = transpose(myModulator(transpose(source_bits)));

guard_bits = (randn(1,10))>0;

tx_frame = [Preamble_1; Preamble_2; mod_symbols'; guard_bits'];

filteredFrme = fltr(tx_frame, s);

%power spectrum
figure()
[Pxx_tx,F_tx] = pwelch(filteredFrme,[],[],[],hw_stp.R,'centered');
plot(F_tx,10*log10(Pxx_tx));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Rx signal')

%99 OBW
figure()
obw(filteredFrme,hw_stp.R,[],99);

%CCDF
figure()
ccdf = comm.CCDF(...
    'AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true);
[ccdfy,ccdfx,avg,peak] = ccdf([filteredFrme filteredFrme]);

pm = powermeter(ComputeCCDF=true);
averagePower = pm(filteredFrme); 
prob = probability(pm,3);
plotCCDF(pm,GaussianReference=true);
 title('CCDF Measurements QPSK')
%constellation diagram
% figure()
% scatterplot(filteredFrme);

%polar diagram
plot(real(filteredFrme),imag(filteredFrme));

%eye diagrams
eyediagram(filteredFrme,16)




%========================================================================

% Filter Setup
s.T ="NON";
s.sps = 8;
s.span = 12;
s.alpha = 0.3;
OvSampRatio = 4;


%gnerate preamble
Preamble_1 = mseq(2,7);

Preamble_2 = mseq(2,7);
%generatiing source bits
source_bits = (randn(1,512))>0;

%qpsk modulation
    myModulator = comm.PSKModulator(...
        'ModulationOrder',  mod_stp.M,...   % here the second input determines how to
        'PhaseOffset',      pi/4, ...      % rotate (phase offset) the constellation
        'BitInput',         true);         % point in the complex plane. If we don't
    % rotate, it provides QPSK costellations
    % on the I and Q axis like [1 j -1 -j].
    mod_symbols = transpose(myModulator(transpose(source_bits)));

guard_bits = (randn(1,10))>0;

tx_frame = [Preamble_1; Preamble_2; mod_symbols'; guard_bits'];

temp            =     repmat(tx_frame,OvSampRatio,1);
filteredFrme    =     reshape(temp,1,length(tx_frame)*OvSampRatio);

%filteredFrme = fltr(tx_frame, s);

%power spectrum
figure()
[Pxx_tx,F_tx] = pwelch(filteredFrme,[],[],[],hw_stp.R,'centered');
plot(F_tx,10*log10(Pxx_tx));
grid on
xlabel('Frequency (Hz)');
ylabel('Normalized Power (dB)')
title('Power Spectrum of the Rx signal')

%99 OBW
figure()
obw(filteredFrme,hw_stp.R,[],99);

%CCDF
figure()
ccdf = comm.CCDF(...
    'AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true);
[ccdfy,ccdfx,avg,peak] = ccdf([filteredFrme filteredFrme]);

pm = powermeter(ComputeCCDF=true);
averagePower = pm(filteredFrme); 
prob = probability(pm,3);
plotCCDF(pm,GaussianReference=true);
 title('CCDF Measurements QPSK')
%constellation diagram
% figure()
% scatterplot(filteredFrme);

%polar diagram
figure()
plot(real(filteredFrme),imag(filteredFrme));

%eye diagrams
eyediagram(filteredFrme,16)



