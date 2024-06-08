% Written by Huseyin Arslan
% 1/05/2008
% Purpose: To teach students a basic digital communication simulation set up.
% Updated: 2/29/2024 by A. Kihero
clear all;
%% DEFINE SIMULATION PARAMETERS
Rsym           = 1000;             % Symbol rate
OvSampRatio    = 16;               % Number of samples per symbol(Duration of the pulse)
Fs             = Rsym*OvSampRatio; % Sampling Frequency
ModOrder       = 4;                % QPSK
MessageLength  = 1000;             % Symbols
SNR            = 10;   
   %% GENERATE MESSAGE BITS (random message)
NumBitsPerSym = log2(ModOrder); % Number of Bits Per Symbol


% Desired signal-to-noise ratio in [dB]
Values = [-10 3 5 10];
numberr = zeros(length(Values),length(Values));
for iter = 1: 1000
BER_1     = zeros(1,length(Values));
    for k =1:length(Values)
      SNR = Values(k);
      nErr = 0;
      source_bits   = randn(1,MessageLength*NumBitsPerSym)>0; % Generate random bits-there are various other ways of doing the same thing
    %%  MODULATION
    if(1) % Our Own Implementation
        % Mapping the generated bits to QPSK symbols:
        % Angle [pi/4 3*pi/4 -3*pi/4 -pi/4] corresponds...
        % to Gray code vector [11 10 00 01], respectively.
        % 1 1 --> (+,+), 1 0 --> (+,-), 0 0 --> (-,-),  0 1 -->(-,+).
        table       = exp(1j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols
        table       = table([0 1 3 2]+1);      % Gray code mapping pattern for QPSK symbols
        inp         = reshape(source_bits,NumBitsPerSym,length(source_bits)/NumBitsPerSym);
        mod_symbols = table([2 1]*inp+1);      % maps transmitted bits into QPSK symbols
    end
    
    temp      = repmat(mod_symbols,OvSampRatio,1);
    tx_signal = reshape(temp,1,length(mod_symbols)*OvSampRatio);
    
    
    len             = length(tx_signal);
    noise_var = 10^(-SNR/10) ;
    noise           = sqrt(noise_var) * (randn(1, len) + 1j*randn(1, len))/sqrt(2);
    
    %% PASSING THE MESSAGE THROUGH AWGN CHANNELsubplo
    
        rx_signal = tx_signal + noise;
        
        Output_Bits = zeros(2,length(mod_symbols));
    
        Down_sampling = reshape(rx_signal,OvSampRatio,length(rx_signal)/OvSampRatio);
    
        Modulated_Symbols  = sum(Down_sampling,1)./size(Down_sampling,1);
    
            for n = 1:length(Modulated_Symbols)
                    %first quadrant error
                        symbol_Error_1 = abs(Modulated_Symbols(n) - table(1));
                    %second quadrant error
                        symbol_Error_2 = abs(Modulated_Symbols(n) - table(2));
                    %third quadrant error
                        symbol_Error_3 = abs(Modulated_Symbols(n) - table(4));
                    %fourth quadrant error
                        symbol_Error_4 = abs(Modulated_Symbols(n) - table(3));
                    
                        minimum_error = min([symbol_Error_1 symbol_Error_2 symbol_Error_3 symbol_Error_4]);
                    
                        if minimum_error == symbol_Error_1 
                               Output_Bits(:,n)         = [0 ; 0];
                        elseif minimum_error == symbol_Error_2
                               Output_Bits(:,n)         = [0 ;1];
                        elseif minimum_error == symbol_Error_3
                               Output_Bits(:,n)         = [1 ;1];
                        elseif minimum_error == symbol_Error_4
                               Output_Bits(:,n)         = [1 ;0];
                        end
            end
            
             Output_Bits = Output_Bits(:)';  %reshape(Output_Bits, 1,2*length(Output_Bits));
            % 
             numberr(k) = biterr(source_bits, Output_Bits);
             BER(k)     = numberr(k)/length(source_bits);
    end
end
 figure
 semilogy(SNR, BER,'k-s','LineWidth',1.5)
 xlabel('SNR (dB)'); ylabel('BER')
 grid on
 box on