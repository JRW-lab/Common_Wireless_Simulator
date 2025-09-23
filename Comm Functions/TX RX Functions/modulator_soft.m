function data_mod_ll = modulator_soft(data_bin_ll, bin_sym_mapping)
% function data_mod_ll = modulator_soft(data_bin_ll, bin_sym_mapping)
%
% soft modulation
% calculate the Log-likelihood of symbols based on the LL information of binary bits
%
% data_bin_ll: 2 x N_bit_per_slot mat, 1st row: ll of 0, 2nd row, ll of 1
% bin_sym_mapping: M_ary*N_bit_per_sym, each row is the binary mapping of 1 sym (in 0's and 1's)
%

[M_ary, N_bit_per_sym] = size(bin_sym_mapping);

if M_ary ~= 2^N_bit_per_sym
	error('M_ary ~= N_bit_per_sym, bin_sym_mapping is in wroing dimension!');
end

N_sym = length(data_bin_ll)/N_bit_per_sym;

data_mod_ll = zeros(M_ary, N_sym);

for m = 1:M_ary
	for k = 1:N_bit_per_sym	;
		current_bit_idx = bin_sym_mapping(m, k)+1; % 1 or 2
		current_ll_vec = data_bin_ll(current_bit_idx, [k:N_bit_per_sym:end]);
		data_mod_ll(m, :) = data_mod_ll(m, :) + current_ll_vec;				
	end
end
