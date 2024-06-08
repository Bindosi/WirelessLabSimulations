clear; clc
OSR = 1;
sys.sps = OSR;
sys.alpha = 0.3;
sys.fs = 1e6;
sys.span = OSR;
NumberOfSymbols = 91;
CPSize  =   8;
fftSize = 64;
load("CDMA12_5KHzSpan2Active.mat");
NumberOfUsers = 4;
m_sq        = mseq(2,6);
prmbl       = [m_sq;m_sq];
prmbl       = prmbl/rms(prmbl);
code0 = [ 1 1 1 1]; 
code1 = [ 1 -1 1 -1]; 
code2 = [ 1 1 -1 -1];
code3 = [ 1 -1 -1 1];

prompt = {'Enter the last two digits of your ID number'};
dlgtitle = 'User Number';
fieldsize = [1 40];
definput = {'eg for Y3230022, put 22'};
opts.Interpreter = 'tex';
studentNumber = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
myNumber = mod(str2num(studentNumber{1}),3)
spreadCodes = [code0;code1;code2;code3];
%% Plotting The received transmission frame
figure
plot(20*log10(abs(Y)));

%% Detecting the edge from all the transmitted frames (Power based edge detection)
    window_size         = sys.sps*16;               %detection window
    absolute_values     = abs(Y).^2;
    vector2matrix       = vec2mat(absolute_values, window_size);
    derivatives         = diff(sum(vector2matrix.'));
    plot(derivatives);
     [~, indx] = max(derivatives);
     plot(derivatives)
     edge_index = indx*window_size;

%% Put 1 in the bracket to check the position of the edge of the selected frame
if (0)
        figure()
        plot(20*log10(abs(Y)));
        hold on
        plot(edge_index , 20*log10(abs(Y)),'r*');
        hold off
end
%% Getting a selected frame from all the received frames of the signal
Symbls_off = 1;
frameStartIndex     = edge_index+Symbls_off;
Framelength         = length(prmbl)*OSR+ NumberOfSymbols*OSR*NumberOfUsers;
frameEndIndex       = frameStartIndex + Framelength + Symbls_off*OSR;
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
    symbols_eq      =   symbols_rx_foc.*conj(est_channel);
    scatterplot(preambles_eq(:));
    scatterplot(symbols_eq(:));
%% despreading and demodulating the received symbols into bit
   userSymbols                = reshape(symbols_eq, NumberOfUsers,length(symbols_eq)/NumberOfUsers)
   deSpreadSymbols            = (userSymbols.' * conj(spreadCodes(myNumber+1,:).'))/NumberOfUsers
   receivedBits               = qamdemod(deSpreadSymbols,2);
   scatterplot(deSpreadSymbols)
   %% converting bits to message
r_bits          = receivedBits(1:NumberOfSymbols)
reshapedBits    = reshape(r_bits,13,7)'
messageBits     = reshape(reshapedBits,91,1)

receivedMessage = char(reshape(bin2dec(reshape(char(messageBits + '0'),7,[]).'),1,[]));

msgbox(receivedMessage,"The received Message")

%% Determining Active Codes
% Number of codes and length of each code
num_codes = size(spreadCodes, 1);
code_length = size(spreadCodes, 2);

% Initialize correlation results
correlation_results = zeros(num_codes, 1);

% Perform cross-correlation between received symbols and each code
for i = 1:num_codes
    % Correlate received symbols with the i-th code
    correlation = xcorr(userSymbols(:), spreadCodes(i, :));
    
    % Since we are interested in the active part, take the maximum absolute correlation value
    correlation_results(i) = max(abs(correlation));
end
stem(1:num_codes, correlation_results); xlim([0 4.5])
xlabel('Code Number')
ylabel('Correlation Value')
title('Correlation Between Spread Codes and Received Symbols')
% Display correlation results
disp('Correlation results for each code:');
disp(correlation_results);

% Determine active codes by finding the two highest correlation values
[~, sorted_indices] = sort(correlation_results, 'descend');
active_codes_indices = sorted_indices(1:2);

% Display active codes
disp('Indices of active codes:');
disp(active_codes_indices);

% Display the active codes
disp('Active codes:');
disp(spreadCodes(active_codes_indices, :));

