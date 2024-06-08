 spread_code=[1 1 1 1 1 1 1 1 ; 1 1 1 1 -1 -1 -1 -1;
 1 1 -1 -1 1 1 -1 -1; 1 1 -1 -1 -1 -1 1 1; 1 -1 1 -1 1 -1 1 -1;
 1 -1 1 -1 -1 1 -1 1; 1 -1 -1 1 1 -1 -1 1; 1 -1 -1 1 -1 1 1 -1];
 % There are 8 orthogonalization codes
 Ncs=352*2; % number of chips in one slot
 Qs=8; % Using only the spreading factor of 8 (orthogonal)
 SNR=20; % desired SNR values in dB
 M=4; %modulation order (QPSK)
 mid=144; %length of midamble
Max_num_users=8; User_all=4;
ord=sort(randperm(8,4));
code_order = spread_code(ord,:);
% the spreading codes that the users use.
% ************** TRANSMITTER ************************
slot=zeros(1,864);
for u=1:User_all % up to 16 users can be multiplexed
% one slot of user data generation
bits(u,:) = randi([0 1],1,(Ncs / Qs)*log2(M));
%******** Modulation with QPSK mapping
symb=qammod(bits(u,:).',M,'InputType','bit',...
'UnitAveragePower',true);
%******** spreading
symb_spread=reshape((symb(:)*code_order(u,:)).',1,Ncs);

chips_midamble(u,:)=(randn(1,144)>0)*2-1;
 %****** Slot formation
 %[DATA(352) MIDAMBLE(144) DATA(352) GUARD(16)]
slot=slot + [symb_spread(1:Ncs/2) chips_midamble(u,:)...
symb_spread(Ncs/2+1:end) zeros(1,16)];
 end
 %******* Channel
 L=2; % Channel impulse response length in chip spaced
 cir=(randn(1,L)+1i*randn(1,L))/sqrt(2)/sqrt(L); % uniform
 % Channel is static and randomly generated over the block.
rx=filter(cir,1,slot); % pass the signal through channel
 noise=(randn(1,length(rx))+1i*randn(1,length(rx)))/sqrt(2);
 sgnl_rx = rx + noise*10^(-SNR/20);% add noise

 for usr_sel = 1:Max_num_users
 for LL=1:L % RAKE RECEIVER (L correlators)
 %******* extract the data from slot (demultiplexing)
 data_r=[sgnl_rx(LL:Ncs/2+LL-1) sgnl_rx(Ncs/2+LL+...
 mid:Ncs/2+LL+mid+Ncs/2-1)];
%******* remove the effect of channel
 sgnl_eq=data_r.*conj(cir(LL));
 % channel tap delays and coefficients are known
 %******* De_spreading
 %sp_code=code_order(usr_sel,:);
 sp_code=spread_code(usr_sel,:);
 x1x=reshape(sgnl_eq,Qs,length(sgnl_eq)/Qs);
 rec_Ds(LL,:) = (x1x.' * conj(sp_code(:)))/Qs;
 end
 s_comb(:,usr_sel) = mean(rec_Ds,1); %combine
 end
 E=10*log10(mean(abs(s_comb).^2)/max(mean(abs(s_comb).^2))); 
 b=bar(E);b(1).BaseValue = -20;