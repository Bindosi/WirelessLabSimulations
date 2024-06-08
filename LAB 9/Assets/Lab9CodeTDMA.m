clear; clc
OSR = 16;
sys.sps = OSR;
sys.alpha = 0.3;
sys.fs = 1e6;
sys.span = OSR;
NumberOfSymbols = 91;
CPSize  =   8;
fftSize = 64;
load("TDMA125KHzSpan.mat");
preambles = mseq(2,6),% figure()

prompt = {'Enter the last two digits of your ID number'};
dlgtitle = 'User Number';
fieldsize = [1 40];
definput = {'eg for Y3230022, put 22'};
opts.Interpreter = 'tex';
studentNumber = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
 myNumber = 3;  % mod(str2num(studentNumber{1}),3)

NumberOfUsers = 4;
m_sq        = mseq(2,6);
prmbl       = m_sq;
prmbl       = prmbl/rms(prmbl);
%% Plotting The received transmission frame
figure
plot(20*log10(abs(Y)));

%% Detecting the edge from all the transmitted frames (Power based edge detection)

window_size         = OSR*16;               % edge detection window
absolute_values     = abs(Y).^2;
vector2matrix       = vec2mat(absolute_values, window_size);
derivatives         = diff(sum(vector2matrix.'));
figure()
plot(derivatives);
[~, indx] = max(derivatives);
edge_index = indx*window_size;

%% Put 1 in the bracket to check the position of the edge of the selected frame
if (0)
        figure()
        plot(20*log10(abs(Y)));
        hold on
        plot(edge_index +199, 20*log10(abs(Y)),'r*');
        hold off
end
%% Getting a selected frame from all the received frames of the signal
Symbls_off =2;
frameStartIndex     = edge_index - Symbls_off*OSR;
Framelength         = length(prmbl)*OSR+ NumberOfSymbols*OSR*NumberOfUsers;
frameEndIndex       = frameStartIndex+Framelength + Symbls_off*OSR;
selected_frame      = Y(frameStartIndex+1:frameEndIndex);

figure
plot(20*log10(abs(selected_frame)));
%% Downsampling the selected frame
    sps = OSR;
    indices = zeros(sps,1);
    for ii = 1:sps
        crr = sum(abs(selected_frame(ii:sps:end)).^2);
        indices(ii) = crr;
    end
    opt                 = find(indices ==  max(indices));
    receivedFrame       = downsample(selected_frame(opt:end),sps);
    preambles_rx        = receivedFrame(1:length(preambles))
    symbols_rx          = receivedFrame(length(preambles)+1:end-1);
    scatterplot(preambles_rx);
    scatterplot(symbols_rx);
%% Finding and correcting frequency offset
    sampx = fftshift(fft(preambles_rx.^2,40*length(preambles)));
    delta_f = (sys.fs)/length(sampx);

    freq_interv=(-(sys.fs))/2:delta_f:((sys.fs)/2)-1;
    figure()
    plot(freq_interv/2, abs(sampx));
    off  = freq_interv(sampx==max(sampx))/2;
    tp = (1:length(preambles_rx))'/sys.fs;
    ts = (1:length(symbols_rx))'/sys.fs;
    % Frequency Offset Correction
    preambles_foc   =   preambles_rx.*exp(-2i*pi*tp*off)
    symbols_rx_foc  =   symbols_rx.*exp(-2i*pi*ts*off)
    scatterplot(preambles_foc)
    scatterplot(symbols_rx_foc)
%% Equalizing the preambles
    est_channel     =   mean(preambles_foc./preambles);
    plot(abs(est_channel))
    preambles_eq    =   preambles_foc./est_channel;
    symbols_eq      =   symbols_rx_foc.*conj(est_channel);
    scatterplot(preambles_eq); 
    scatterplot(symbols_eq); 
%% demodulating the received symbols into bit
   receivedBits            = qamdemod(symbols_rx,2);
   %% converting bits to message
startingIndex   = NumberOfSymbols*myNumber+1;
endingIndex     = startingIndex+NumberOfSymbols-1;

r_bits          = receivedBits(startingIndex:endingIndex);
reshapedBits    = reshape(r_bits,13,7)';
messageBits     = reshape(reshapedBits,91,1);

receivedMessage = char(reshape(bin2dec(reshape(char(messageBits + '0'),7,[]).'),1,[]));

msgbox(receivedMessage,"The received Message")
