%Effect of the windowing operation on spectrum of an OFDM signal
%Please rerun the script with different roll-off factors
N = 256;
CP = 1/4;
G = CP*N;

%RC windowing CP
CP2 = 1/16; % the ratio of the additional guard period
G2 = CP2*N;

num_OFDM_symbols=100;

frame=[]; frame_RC=[]; x0 = zeros(1, N+G);
%------Raised-cosine (comparison)-------------
RRC = 0.5 + 0.5.*cos(pi + pi.*(0:G2-1)./G2);
figure()
plot(RRC)
LRC = 0.5 + 0.5.*cos(pi.*(0:G2-1)./G2);
figure()
plot(LRC)
%---------------------------------------------

for ii=1:num_OFDM_symbols
%symbols=2*(randn(1,N)>0)-1;
symbols= qammod( floor(rand(1,N)*4) ,4);
used_subc_indices=[floor(N/4):1:N/2-1 N/2+1:1:N-floor(N/4)];
to_ifft=zeros(1,N);
to_ifft(used_subc_indices)=symbols(used_subc_indices);
time=ifft(to_ifft);
time_cp=[time(end-G+1:end) time];
frame=[frame time_cp];

%------Raised-cosine (comparison)-------------
Ccp = LRC.*x0(G+[1:G2]) + RRC.*time_cp(end-G-G2+[1:G2]);

frame_RC = [frame_RC, Ccp, time_cp];

x0 = time_cp;

end
figure(1)
clf
[Pxx,F]=pwelch(frame,[],[],2048,1,'twosided');
plot(F,10*log10(Pxx))
hold on
[Pxx,F]=pwelch(frame_RC,[],[],2048,1,'twosided');
plot(F,10*log10(Pxx),'r--');