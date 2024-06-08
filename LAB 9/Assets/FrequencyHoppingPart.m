% Frequency Hopping spread spectrum
clear; clc
load("FreqHoping125KhzSpan.mat");

%% Plotting The received transmission frame
figure
plot(20*log10(abs(Y)));

d = seconds(6.25e-6);
win = hamming(100,"periodic");

stft(Y,d,Window=win,OverlapLength=98,FFTLength=128);