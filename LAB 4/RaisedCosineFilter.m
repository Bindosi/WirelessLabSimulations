

index = 0:2;
phase = -i*(2*pi*index/3);

x = sum(exp(phase));

m = 4;
data_symbols = randi([0 m-1], 512, 1)>0;
mod_symbols = qammod(data_symbols, m);

guard_symbols = randi();

tx_frame = [preamble_1 preamble_2 mod_symbols, guard_symbols];
alpha = [0.0 0.3 0.5 0.7 1.0];

filter = rcosdesign(tx_frame, flt.stp.span, flt.stp.sps,"normal");