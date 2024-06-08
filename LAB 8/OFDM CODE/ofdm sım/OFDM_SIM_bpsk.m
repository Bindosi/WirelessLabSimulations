clc
clear all
close all


%% parameters
N=64;% number of subcar
data=2*(randi([0 1],N,1))-1;
st=ifft(data);
figure
plot(abs(st),'linewidth',3)
lcp=4;
st_cp=[st(end-lcp:end);st];% add cp

%% channel
taps=3;
h=randn(taps,1)+1i*randn(taps,1);
y=conv(st_cp,h);
y=awgn(y,5,'measured');
%% receiver 
y=y(lcp+1:N+lcp);
yf=fft(y);

H=fft(h,N);
figure
plot(abs(H),'linewidth',3)
H_est=fft(ifft(yf(1:8:end)./data(1:8:end)),N);
figure
plot(abs(H_est),'LineWidth',3)
hold on
plot(abs(H),'*')
X_dem=yf./H_est;
%X_dem=2*(X_dem>0)-1;
X_dem=real(X_dem)
for i=1:N
    if X_dem(i)<0
        X_dem(i)=-1
    else
        X_dem(i)=1
    end
end
       
scatterplot(X_dem)
scatterplot(data)
ber=sum(data~=X_dem)/N







