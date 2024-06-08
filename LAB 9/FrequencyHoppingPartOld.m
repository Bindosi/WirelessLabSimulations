% Frequency Hopping spread spectrum
clear; clc
OSR = 1;
sys.fs = 1e6;
NumberOfSymbols = 91;
load("FreqHoping125KhzSpan.mat");
myNumber = mod(21,3)
NumberOfUsers = 4;


d = seconds(1e-6);
win = hamming(100,"periodic");

stft(Y,d,Window=win,OverlapLength=98,FFTLength=128);

sps = 16;           % Samples per symbol
channelBW = 2e6;    % Channel spacing (Hz) as per standard
symbolRate = 1e6;
x = [];
winL = 91; noverlap = 50; nfft = 64;
sampleRate = symbolRate*sps*10;

[s,f,t]=spectrogram(Y,winL,noverlap,nfft,sampleRate,'yaxis');
imagesc(t/1e-6,f/1e6,10*log10(abs(s)))
set(gca,'YDir','normal'); cc = colorbar;
ylabel(cc,'Power/frequency (dB/Hz)');
hColourbar.Label.Position(1) = 3;
