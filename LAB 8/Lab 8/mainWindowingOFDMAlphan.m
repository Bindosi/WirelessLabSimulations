%==========================================================================
% Alphan Sahin, 2015
% Windowing OFDM
%==========================================================================
clc;
clear all;
close all;

% Main Parameters
W = 4;
N = 512;	% Number of subcarriers for DFT
G = 64;    	% Length of CP
Nblk = 300; % Number of OFDM symbols considered in the simulations
M = 16;     % M-QAM



% OFDM Settings
F = dftmtx(N)/sqrt(N);
A = [zeros(G,N-G) eye(G); eye(N)]; % CP addition matrix

%% Simulation
% Transmitter
indices = 1:100;

%QAM
m = 16;
data = randi([0,1], numel(indices)*log2(M),Nblk);
mod_symbols16 = qammod(data, M ,UnitAveragePower=true);

d = zeros(N,Nblk);
d(indices,:) = data(indices,:);
xOFDM = A*F'*d;

w1 = 1/2*(cos(linspace(-pi,0,W+2))+1);
plot(w1)

w1(1) = [];
w1(end) = [];

W1 = eye(N+G);
W1(1:W,1:W) = diag(w1);

w2 = 1/2*(cos(linspace(0,pi,W+2))+1);
plot(w2)

w2(1) = [];
w2(end) = [];
W2 = zeros(N+G);

W2(1:W,(G+(1:W))) = diag(w2);

xWOFDM = W1*xOFDM  + W2* xOFDM;

if (1) % OOB
    xA = xOFDM(:);
    Nwin = 1024;
    Noverlap= floor(Nwin*3.2/4);
    [pX,f] = pwelch(xA(:), Nwin, Noverlap,[],N);
    pX = pX * N;
    h100 = figure(100);
    plot(f-N/2,fftshift(10*log10(pX)) , 'k-', 'displayname', 'OFDM' )
    hold on

    xA = xWOFDM(:).';
    Nwin = 1024;
    Noverlap= floor(Nwin*3.2/4);
    [pX,f] = pwelch(xA(:), Nwin, Noverlap,[],N);
    pX = pX * N;
    h100 = figure(100);
    plot(f-N/2,fftshift(10*log10(pX)) , 'k--', 'displayname', 'Windowing & OFDM' )
    hold on
    
    
    grid on
    xlabel('Subcarrier Index')
    ylabel('Power spectral density [dBm/Hz]')
    ylim([-50 5])
    xlim([-N/2,N/2-1])
    drawnow
    legend('show')
    legend('Location','Best')
end

