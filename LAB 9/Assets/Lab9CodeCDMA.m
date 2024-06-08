clear; clc
OSR = 1;
sys.sps = OSR;
sys.alpha = 0.3;
sys.fs = 1e6;
sys.span = OSR;
NumberOfSymbols = 91;
CPSize  =   8;
fftSize = 64;
load("CDMA12_5KHzSpan.mat");

prompt              = {'Enter the last two digits of your ID number'};
dlgtitle            = 'User Number';
fieldsize           = [1 40];
definput            = {'eg for Y3230022, put 22'};
opts.Interpreter    = 'tex';
studentNumber       = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
myNumber            = mod(str2num(studentNumber{1}),3)
NumberOfUsers       = 4;
m_sq                = mseq(2,6);
prmbl               = [m_sq;m_sq];
prmbl               = prmbl/rms(prmbl);
code0 = [ 1 1 1 1]; 
code1 = [ 1 -1 1 -1]; 
code2 = [ 1 1 -1 -1];
code3 = [ 1 -1 -1 1];

spreadCodes = [code0;code1;code2;code3];
%% Plotting The received transmission frame
figure
plot(20*log10(abs(Y)));

%% Detecting the edge from all the transmitted frames (Correlation based edge detection)
PreamblesCreated    =   [mseq(2,6); mseq(2,6)];
correlation         =   xcorr(PreamblesCreated,Y);
[~,indx1]           =   maxk(correlation,5);
sortedIndexes       =   sort(indx1);

figure()
plot(abs(correlation))

edge_index      = min(sortedIndexes);
%% Put 1 in the bracket to check the position of the edge of the selected frame
if (0)
        figure()
        plot(20*log10(abs(Y)));
        hold on
        plot(edge_index+53 , 20*log10(abs(Y)),'r*');
        hold off
end
%% Getting a selected frame from all the received frames of the signal
Symbls_off = 53;
frameStartIndex     = edge_index+Symbls_off;
Framelength         = length(prmbl)*OSR+ NumberOfSymbols*OSR*NumberOfUsers;
frameEndIndex       = frameStartIndex +1+ Framelength;
selected_frame      = Y(frameStartIndex+1:frameEndIndex);

%% Downsampling the selected frame
    preambles_rx        = selected_frame(1:length(prmbl))
    symbols_rx          = selected_frame(length(prmbl)+1:end-1);

    scatterplot(preambles_rx);
    scatterplot(symbols_rx);
    figure
    plot(20*log10(abs(selected_frame)));
    hold on
    plot(20*log10(abs(preambles_rx)),"red");
%% Finding and correcting frequency offset
    phase=angle(sum(prmbl.*conj(preambles_rx)));
    preambles_foc   =   preambles_rx.*exp(i*phase);
    symbols_rx_foc  =   symbols_rx.*exp(i*phase);
    scatterplot(preambles_foc);
    scatterplot(symbols_rx_foc);
%% Equalizing the preambles
    est_channel     =   mean(preambles_foc./prmbl);
    preambles_eq    =   preambles_foc./est_channel;
    symbols_eq      =   symbols_rx.*conj(est_channel);
    scatterplot(preambles_eq(:));
    scatterplot(symbols_eq(:));
%% despreading and demodulating the received symbols into bit
   userSymbols                = reshape(symbols_eq, NumberOfUsers,length(symbols_eq)/NumberOfUsers)
   deSpreadSymbols            = (userSymbols.' * conj(spreadCodes(myNumber+1,:).'))/NumberOfUsers
   receivedBits               = qamdemod(deSpreadSymbols,2);
  
   %% converting bits to message
r_bits          = receivedBits(1:NumberOfSymbols)
reshapedBits    = reshape(r_bits,13,7)'
messageBits     = reshape(reshapedBits,NumberOfSymbols,1)

receivedMessage = char(reshape(bin2dec(reshape(char(messageBits + '0'),7,[]).'),1,[]));

msgbox(receivedMessage,"The received Message")

sum1 = abs(sum((userSymbols.' * conj(code0(:)))))
sum2 = abs(sum((userSymbols.' * conj(code1(:)))))
sum3 = abs(sum((userSymbols.' * conj(code2(:)))))
sum4 = abs(sum((userSymbols.' * conj(code3(:)))))
