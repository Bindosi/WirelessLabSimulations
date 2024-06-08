clc
clear all
close all
%% parameters
N_FFT       =   64;% number of subcar
M           =   4;
data        =   (randi([0 M-1],N_FFT,1));
mod_data    =   qammod(data,M)

% scatterplot(mod_data)
st          =   ifft(mod_data);

figure
plot(abs(st),'linewidth',3)
lcp=4;
st_cp=[st(end-lcp:end);st];% add cp

%% channel
iter=4000;
for i=1:iter
taps=3;
h=randn(taps,1)+sqrt(-1)*randn(taps,1);
SNR=0:5:30;
for s=1:length(SNR)
y=conv(st_cp,h);
y=awgn(y,SNR(s),'measured');
%% receiver 
y=y(lcp+1:N_FFT+lcp);
yf=fft(y);
H=fft(h,N_FFT);
% figure
% plot(abs(H),'linewidth',3)
H_est=fft(ifft(yf(1:4:end)./mod_data(1:4:end)),N_FFT);
% figure
% plot(abs(H_est),'LineWidth',3)
% hold on
% plot(abs(H),'*')
X_dem=yf./H_est;% equalize
% scatterplot(X_dem)
X_dem=qamdemod(X_dem,M);
%% please convert this back to binary (both data and X_dem) then do BER
ber(s,i)=sum(data~=X_dem)/N_FFT;
end
end
avg_ber=mean(ber.');
figure
semilogy(SNR,avg_ber,'LineWidth',3)
grid on
















