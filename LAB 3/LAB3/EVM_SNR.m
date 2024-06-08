
data_symbols16 = randi([0 m-1], 512, 1);
mod_symbols16 = qammod(data_symbols16, m ,UnitAveragePower=true);


evm = comm.EVM;
rmsEVM4 = evm(temp4qm(end-511:end),symbs_rx4qm/rms(symbs_rx4qm));

scatterplot(temp4qm(end-511:end))
hold on
scatterplot(symbs_rx4qm/rms(symbs_rx4qm))

% SNR4 = 10*log10(signalp/signalp_noise)
av4=mean(abs(symbs_rx4qm).^0.5);
rmsEVM4 =rmsEVM4/av4; 
SNR_EVM4 = -20*log10(rmsEVM4/100);

evm = comm.EVM;
rmsEVM16 = evm(mod_symbols16,symbs_rx16qm/rms(symbs_rx16qm));



















scatterplot(temp16qm(end-511:end))
hold on
scatterplot(symbs_rx16qm/rms(symbs_rx16qm))

% SNR16 = 10*log10(signalp/signalp_noise)
av16=mean(abs(symbs_rx16qm).^0.5);
rmsEVM16 =rmsEVM16/av16; 
SNR_EVM = -20*log10(rmsEVM16/100);

evm = comm.EVM;
rmsEVM64 = evm(temp64qm(end-511:end),symbs_rx64qm/rms(symbs_rx64qm));

scatterplot(temp64qm(end-511:end))
hold on
scatterplot(symbs_rx64qm/rms(symbs_rx64qm))

% SNR = 10*log10(signalp/signalp_noise)
av64=mean(abs(symbs_rx64qm).^0.5);
rmsEVM64 =rmsEVM64/av64; 
SNR_EVM64 = -20*log10(rmsEVM64/100);

% Q3 contellation and Eye Diagrams

