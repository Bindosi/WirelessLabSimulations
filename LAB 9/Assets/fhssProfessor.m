numPackets = 10; % Number of packets to generate
sps = 16; % Samples per symbol
 messageLen = 2000; % Length of message in bits
phyMode = 'LE1M'; % Select mode
channelBW = 2e6; % Channel spacing (Hz) as per standard
symbolRate = 1e6;
 x = [];
 winL = 1000; noverlap = 500; nfft = 8192;
sampleRate = symbolRate*sps*10;
 % Loop over the number of packets, generate a BLE waveform
 rng default;
 for packetIdx = 1:numPackets
 message = randi([0 1],messageLen,1); % Message bits
 chanIndex = randi([0 39],1,1); % Channel index
 if(chanIndex >=37)
% Access address for periodic advertising channels
accessAdd = [0 1 1 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 0 ...
1 0 0 0 1 0 1 1 1 0 0 0 1]';
 else
 % Random access address for data channels
 accessAdd = [0 0 0 0 0 0 0 1 0 0 1 0 0 ...
 0 1 1 0 1 0 0 0 1 0 1 0 1 1 0 0 1 1 1]';
 end
 waveform = bleWaveformGenerator(message,'Mode',phyMode,...
 'SamplesPerSymbol',sps,...
 'ChannelIndex',chanIndex,...
 'AccessAddress',accessAdd);
 frequencyOffset = channelBW*chanIndex;
 waveform = waveform.* exp(1i*2*pi*frequencyOffset*...
 (1:length(waveform)).'/sampleRate);
 x = [x; waveform];
 end
 [s,f,t]=spectrogram(x,winL,noverlap,nfft,sampleRate,'yaxis');
 imagesc(t/1e-6,f/1e6,10*log10(abs(s))); ylim([0, 100])
 set(gca,'YDir','normal'); cc = colorbar;
 ylabel(cc,'Power/frequency (dB/Hz)');
 hColourbar.Label.Position(1) = 3