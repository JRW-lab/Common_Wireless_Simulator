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

% encoder'*encoder does not depend on k (encoder is fixed for this whole
% call) - computed once here instead of recomputed via a full matrix
% product on every pass of the k-loop below.
EtE = encoder' * encoder;

decoded_soft = zeros(2, N_tx_sym);
decoded_hard = zeros(N_tx_sym, 1);
%decoded_hard(N_tx_sym) = modulator(demodulator(decoded_soft(N_tx_sym), M_ary), M_ary);

% ---------------------------------------------------------------------
% Rank-1-update precomputation. For every k, Rss only differs from
% diag(data_var) at position k (set to 1 there), so both
% Ree(k) = inv_Rss(k) + EtE/noise_var and Rnn+encoder*Rss(k)*encoder'
% are each a COMMON, k-independent base matrix plus a single rank-1
% correction at position/column k. This lets the Cholesky factor (via
% cholupdate, O(N^2)) and the MMSE matrix's inverse (via
% Sherman-Morrison, O(N^2)) be derived per k from one base
% decomposition/inverse computed once, instead of a full O(N^3)
% chol()/matrix-solve recomputed from scratch for every one of the
% N_tx_sym values of k. Verified equivalent via
% scratchpad/validate_rank1_math.m before this rewrite.
data_var_floored = data_var;
data_var_floored(data_var_floored < 1e-6) = 1e-6;

c_vec = 1 - 1./data_var_floored;      % Ree(k) = Ree_base + c_vec(k)*e_k*e_k'
d_vec = 1 - data_var_floored;         % M(k)   = M_base   + d_vec(k)*v_k*v_k', v_k = encoder(:,k)

Ree_base = diag(1./data_var_floored) + EtE/noise_var;
Ree_base(logical(eye(size(Ree_base)))) = abs(diag(Ree_base));
R_base = chol(Ree_base);

M_base = Rnn + (encoder .* data_var_floored) * encoder';
Minv_base = inv(M_base);
W = Minv_base * encoder;                      % W(:,k) = Minv_base*encoder(:,k)
alpha_vec = real(sum(conj(encoder).*W,1)).';   % alpha_vec(k) = encoder(:,k)'*W(:,k)
% ---------------------------------------------------------------------

d0_mat = zeros(N_bit_per_sym, N_half_state);
d1_mat = zeros(N_bit_per_sym, N_half_state);
for k = N_tx_sym:-1:1

    mean_vec = data_mean;
    mean_vec(k) = 0;

    % Cholesky row k via rank-1 update/downdate of the base factor,
    % falling back to a direct from-scratch decomposition of Ree(k) in
    % the rare case the update/downdate is numerically unsafe.
    c_k = c_vec(k);
    if c_k == 0
        R_k = R_base;
    else
        x = zeros(N_tx_sym,1);
        x(k) = sqrt(abs(c_k));
        try
            if c_k > 0
                R_k = cholupdate(R_base, x, '+');
            else
                R_k = cholupdate(R_base, x, '-');
            end
            if ~isfinite(R_k(k,k)) || ~isreal(R_k(k,k)) || R_k(k,k) <= 0
                error('bdfe:choldowndate','unsafe update');
            end
        catch
            Ree_k = Ree_base;
            Ree_k(k,k) = Ree_k(k,k) + c_k;
            R_k = chol(Ree_k);
        end
    end
    r_kk = R_k(k,k);
    U_row = R_k(k,:) / r_kk;   % U(k,:): unit-diagonal normalized row
    Dkk = r_kk^2;              % D(k,k)

    % MMSE row via Sherman-Morrison, computed as a chain of
    % matrix-vector products so a full N_tx_sym x N_rx_sym MMSE matrix
    % is never formed (U(k,k) is always 1 by construction, so the
    % Feedback_mat(k,k) multiplier from the original code is omitted -
    % multiplying by 1 exactly, not an approximation).
    var_row = data_var_floored;
    var_row(k) = 1;
    p_k = (U_row .* var_row) * encoder';

    d_k = d_vec(k);
    w_k = W(:,k);
    denom = 1 + d_k*alpha_vec(k);
    Feedforward_row = (Minv_base*p_k' - (d_k/denom)*w_k*(w_k'*p_k'))';

    data_ff(k) = Feedforward_row*(data_rcv.'-encoder*mean_vec.');
    temp_value = data_ff(k)-U_row(k+1:end)*(decoded_hard(k+1:end)-data_mean(k+1:end).');
    dec_vec = -abs(temp_value-syms.').^2*Dkk;
    % find the LL of '0' and '1'
    for m = 1:N_bit_per_sym
        d0_mat(m, :) = -abs(temp_value-group0(m, :)).^2*Dkk;
        d1_mat(m, :) = -abs(temp_value-group1(m, :)).^2*Dkk;
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
