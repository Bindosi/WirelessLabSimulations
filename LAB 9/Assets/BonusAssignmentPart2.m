%% Code For Bonus Assignment
% Registration Number D3230014
% Turkish ID = 98674047230

clear;
clc;
c1 = 4; c2 = 2; 

% Since Taking c1/10 will make the term 0.9 + c1/10 greater than one and
% thus making the statement 0.9+c1/10 < |H(e^jw)| <= 1 for 0<w<0.2pi
% invalid, we take it as 0.9/100

HejwUp      = (0.9+c1/100)^2;
HejwLow     = (0.1+c2/100)^2;

N = 2;

Wp = 2*tan(0.1*pi); 
Ws = 2*tan(0.25/2*pi);

Wc = 2*tan(Wp+Ws)/2;
tol         = 10000;
while tol>1e-6
    b       = 1/HejwUp;
    c       = b-1;
    d       = log(c)/(2*N);
    WCtemp  = 0.2*pi/exp(d);
    b       = 1/HejwLow;
    c       = b-1;
    N       = log(c)/log(0.3*pi/WCtemp)/2;
    tol     = abs(Wc-WCtemp);
    Wc      = WCtemp;
end

filterOrder = 2*ceil(N);

% We then determine the value of Wc for the obtained value of Nc
WcforNc         = 0.3*pi/nthroot((1/HejwLow - 1),filterOrder);


% Here we determine the poles of C(s) = H(s) * H(-s)
a               = 1;
b               = 1/(WcforNc^filterOrder);
CsVector        = [b, zeros(1, filterOrder-1), a];
CsRoots         = roots(CsVector);

figure()
zplane([],CsRoots);

% We then pick the roots on the left of the s plane for stability
stableRoots =  CsRoots(real(CsRoots)<0);

% Checking pole placement
figure();
zplane([],stableRoots);

% We then determine the value of K for H(0) = 1 and S = 0
K               = 1;
for root        = 1:length(stableRoots)
K               = -K * (stableRoots(root));
end

numeS                   = K;
denomeS                 = poly(stableRoots); 
[sCoeff, sPoles, k]     = residue(numeS,denomeS);

% To get the the poles of Z we use (1-exp(ps)Z^-1) where ps =  poles in S domain
zPoles                  = exp(sPoles);
[num, den]              = residue(sCoeff, zPoles, k);

[H,w]=freqz(num,den);
figure
title('Magnitude Hejw against frequency in radians');
xlabel('Frequency w');
ylabel('Magnitude Hejw w');
plot(w/pi,abs(H));



%% Checking the filter with butterworth matlab function
clear all;
close all;
 Rp = 0.25;
 Rs = 30;
 Wp = [0.2];
 Ws = [0.25];

 FilterType = 'Butterworth Filter';

 [N, Wpe] = buttord(Wp, Ws, Rp, Rs);
 [b,a] = butter(N,Wpe);
 [H,w] = freqz(b,a);
 %[Gd] = grpdelay()

 figure
plot(w/pi,20*log10(abs(H)))
title('Magnitude Hejw against frequency in radians');
xlabel('Frequency w');
ylabel('Magnitude Hejw w');



