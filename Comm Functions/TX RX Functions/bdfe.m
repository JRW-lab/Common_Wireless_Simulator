function [decoded_hard, ll_bit_extrinsic] = bdfe(data_rcv, ll_bit_apriori, encoder, noise_var, group0, group1, bin_mat, syms)
% function decoded = bdfe(data_rcv, llr_bit, encoder, noise_var, N_bit_per_sym, group0, group1, bin_mat, Feedforward_mat, Feedback_mat)
%   Block Decision Feedback Equalizer
%	Stamoulis TCOM 2001
%
%	softoutput: log likelihood
%		    ln(Pr(x = 1)); ln(Pr(x = 0))
%
% 	model:
%	data_rcv = encoder*data + noise
%
%	input:
%	ll_bit_apriori: a priori Log-likelihood of the data from the prevous iteration
%	noise_var: variance of AWGN
%	group0 and group1: please refer to modulator.m
%	bin_mat: the matrix contains the unmodulated bits. Please refer to the main program
%	syms: the modulated symbols.
%


M_ary = length(syms);
N_bit_per_sym = log2(M_ary);

[temp_row, temp_col] = size(bin_mat);

if or(M_ary ~= temp_row, N_bit_per_sym ~= temp_col)
    error('bin mat and syms mismatch');
end

[N_bit_per_sym, N_half_state] = size(group0);
if N_half_state ~= 2^N_bit_per_sym/2
    error('size of group0 is wrong');
end

% conver bit a priori information to symbol a priori information
% (convert bit LL to sym LL)
ll_sym_apriori = modulator_soft(ll_bit_apriori, bin_mat);


[N_rx_sym, N_tx_sym] = size(encoder);

% normalize the a priori information
log_normalize_coeff_vec = zeros(1,N_tx_sym);
for m = 1:N_tx_sym
    log_normalize_coeff_vec(m) = 0-logsum(ll_sym_apriori(:, m).');
end
log_sym_apriori = ll_sym_apriori + ones(M_ary, 1)*log_normalize_coeff_vec;
p_sym_apriori = exp(log_sym_apriori);
% p_sym_apriori = p_sym_apriori./(ones(M_ary, 1)*sum(p_sym_apriori, 1));
p_sym_apriori = p_sym_apriori ./ sum(p_sym_apriori,1);


% p_bit_apriori = exp(ll_bit_apriori);
% p_bit_apriori = p_bit_apriori./(ones(2, 1)*sum(p_bit_apriori, 1));

% group0_idx = zeros(N_bit_per_sym,N_half_state);
% group1_idx = zeros(N_bit_per_sym,N_half_state);
% for m = 1:N_bit_per_sym
% 	for k = 1:N_half_state
% 		group0_idx(m, k) = (group0(m, k) == syms)*(1:length(syms)).';
% 		group1_idx(m, k) = (group1(m, k) == syms)*(1:length(syms)).';
% 	end
% end
[~, group0_idx] = ismember(group0, syms);
[~, group1_idx] = ismember(group1, syms);

% mod_level = 2^N_bit_per_sym;

% Es = 1;

% calculate the mean and variance of the symbols
% mean
data_mean = syms*p_sym_apriori;
data_var = abs(syms).^2*p_sym_apriori-abs(data_mean).^2;


% correlation matrix of the tx symbols
% Rss = diag(data_var);
data_var(data_var == 0) = 1e-5;
% inv_Rss = diag(1./data_var);

% correlation matrix of AWGN
Rnn = noise_var*eye(N_rx_sym);

decoded_soft = zeros(2, N_tx_sym);
decoded_hard = zeros(N_tx_sym, 1);
%decoded_hard(N_tx_sym) = modulator(demodulator(decoded_soft(N_tx_sym), M_ary), M_ary);


d0_mat = zeros(N_bit_per_sym, N_half_state);
d1_mat = zeros(N_bit_per_sym, N_half_state);
for k = N_tx_sym:-1:1

    % start of matrix formulation
    % correlation matrix of the tx symbols
    var_vec = data_var;
    var_vec(k) = 1;
    mean_vec = data_mean;
    mean_vec(k) = 0;
    Rss = diag(var_vec);
    % debug
    var_vec(var_vec < 1e-6) = 1e-6;
    inv_Rss = diag(1./var_vec);


    % MMSE matrix
    % if rcond(Rnn+encoder*Rss*encoder') > 1e-2
    % 	MMSE_mat = Rss*encoder'/(Rnn+encoder*Rss*encoder');
    % else
    % 	MMSE_mat = Rss*encoder'*pinv(Rnn+encoder*Rss*encoder');
    % end
    MMSE_mat = (Rss * encoder') * ((Rnn + encoder*Rss*encoder') \ eye(N_rx_sym));

    % formulate the feedbackword matrix
    % U'*D*U = Ree, where U is upper triangular with unit diagonal
    Ree = inv_Rss + 1/noise_var*eye(N_tx_sym)*(encoder'*encoder);
    Ree(logical(eye(size(Ree)))) = abs(diag(Ree));
    U_temp = chol(Ree);
    D = diag(diag(U_temp).^2);
    norm_mat = diag(U_temp)*ones(1, N_tx_sym);
    U = U_temp./norm_mat;

    % Feedback_mat = U - eye(N_tx_sym);
    Feedback_mat = U;
    Feedforward_mat = U*MMSE_mat;
    % end of matrix formulation


    data_ff(k) = Feedforward_mat(k, :)*(data_rcv.'-encoder*mean_vec.');
    temp_value = data_ff(k)-Feedback_mat(k, k+1:end)*(decoded_hard(k+1:end)-data_mean(k+1:end).');
    dec_vec = -abs(temp_value-Feedback_mat(k, k)*syms.').^2*D(k, k);
    % find the LL of '0' and '1'
    for m = 1:N_bit_per_sym
        d0_mat(m, :) = -abs(temp_value-Feedback_mat(k, k)*group0(m, :)).^2*D(k, k);
        d1_mat(m, :) = -abs(temp_value-Feedback_mat(k, k)*group1(m, :)).^2*D(k, k);
    end

    for m = 1:N_bit_per_sym
    	ll_app0 = d0_mat(m, :) + ll_sym_apriori(group0_idx(m, :), k).';
    	ll_app1 = d1_mat(m, :) + ll_sym_apriori(group1_idx(m, :), k).';
    	decoded_soft(1, (k-1)*N_bit_per_sym+m) = logsum(ll_app0);
    	decoded_soft(2, (k-1)*N_bit_per_sym+m) = logsum(ll_app1);
    end



    dec_vec = dec_vec+ll_sym_apriori(:, k);
    % normalize
    dec_vec = dec_vec - logsum(dec_vec.');
    decoded_hard(k) = syms*exp(dec_vec);

end
decoded_hard = decoded_hard.';


ll_bit_extrinsic = decoded_soft - ll_bit_apriori;
% debug
%ll_bit_extrinsic = decoded_soft;

% normalize
llr_bit_extrinsic = ll_bit_extrinsic(2, :)-ll_bit_extrinsic(1, :);
exp_llr = exp(llr_bit_extrinsic);
isinf_flag = isinf(exp_llr);
common_llr_vec(~isinf_flag) = log(1 + exp_llr(~isinf_flag));
common_llr_vec(isinf_flag) = llr_bit_extrinsic(isinf_flag);

ll_bit_extrinsic(1, :) = -common_llr_vec;
ll_bit_extrinsic(2, :) = llr_bit_extrinsic-common_llr_vec;
