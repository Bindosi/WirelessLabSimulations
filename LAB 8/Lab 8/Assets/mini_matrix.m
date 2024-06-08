clear
clc

cp          = 2;
guard       = 1;
nfft_size   = 4;
a           = [2 2 2 2;4 4 4 4;6 6 6 6;8 8 8 8;]

w2          = [0 0 0 0;0 0 0 0.5;0 0 0 0;0 0 0 0;]

w1          = [0.5 0 0 0;0 0.5 0 0;0 0 1 0;0 0 0 1;]


cp_addition_matrix = [zeros(nfft_size-cp) eye(cp); eye(nfft_size)]
%cp_addition_matrix  = [0 0 1 0;
%                       0 0 0 1;
%                       1 0 0 0;
%                       0 1 0 0;
%                       0 0 1 0;
%                       0 0 0 1;]

a_cp    = cp_addition_matrix*a; 

add_guard_low = a_cp(length(a_cp)-cp-guard+1:length(a_cp)-cp,:) ;
add_guard_high = a_cp(cp+1,:);

matrix_witLowerAndHigherGuard = [add_guard_low; a_cp;add_guard_high]

windowed = w1*matrix_witLowerAndHigherGuard + w2*matrix_witLowerAndHigherGuard


fftmatrix = [1;2;3;4;5;1;2;3;4;5;1;2;3;4;5;1;2;3;4;5;1;2;3;4;5;1;2;3;4;5]

reshapedFFT = reshape(fftmatrix, 5,6)


pilotBits             = repmat([1 0],1,8);
pilotBitsModulated    = pilotBits*2-1;
pilotBitsReshaped     = repmat(pilotBitsModulated(:),1,25);




