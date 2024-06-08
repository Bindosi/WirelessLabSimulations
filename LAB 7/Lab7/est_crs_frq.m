function off= est_crs_frq(smpls,M,Pr,fs)
M=M-Pr;

smpls_pwr=smpls.*M;
smpls_fft=(abs(fftshift(fft(smpls_pwr))));
df=fs/length(smpls_fft);
df=fs/length(smpls_fft);
freq=-fs/2:df:fs/2-1;
off=freq(smpls_fft==max(smpls_fft))/M;
end