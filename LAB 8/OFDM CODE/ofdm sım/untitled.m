clc
clear all
close all

subcarrier_spacing  =   15e3;
active_subcarriers  =   48;
fft_size            =   64;
number_of_symbols   =   25;
oversampling_rate   =   8;
filter_type         = 'rectangular';
tau_max             =   3.2*10^-5;
sampling_rate       =  2*subcarrier_spacing;
taps                =   ceil(subcarrier_spacing*tau_max);
cp_size             =   taps;

