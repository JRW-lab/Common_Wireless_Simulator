function [decoded_hard, decoded_soft] = bdfe_OFDM(data_rcv, encoder, noise_variance, S)
% function decoded = bdfe_OFDM(data_rcv, encoder, noise_variance, N_bit_per_sym)
%   Block Decision Feedback Equalizer (OFDM, hard-decision only)
%	Stamoulis TCOM 2001
%
% Renamed from bdfe.m: this filename collided with the unrelated 8-arg
% turbo/soft-decision Comm Functions/TX RX Functions/bdfe.m (used by
% sim_fun_OTFS.m). MATLAB path order made TX RX Functions/bdfe.m always
% shadow this one, so sim_fun_OFDM_v2.m's "DD-BDFE" receiver was silently
% calling the wrong function and crashing with "Not enough input
% arguments."

mod_level = length(S);

[N_rx_sym, N_tx_sym] = size(encoder);

% correlation matrix of the tx symbols
Es = 1;
% Rss = Es*eye(N_tx_sym);

% correlation matrix of AWGN
Rnn = noise_variance*eye(N_rx_sym);

% MMSE matrix
MMSE_mat = Es*encoder'*inv(Rnn+Es*encoder*encoder');

% formulate the feedbackword matrix
% U'*D*U = Ree, where U is upper triangular with unit diagonal
Ree = 1/Es*eye(N_tx_sym)+1/noise_variance*eye(N_tx_sym)*encoder'*encoder;
Ree(logical(eye(size(Ree)))) = abs(diag(Ree));
U_temp = chol(Ree);
D = diag(diag(U_temp).^2);
norm_mat = diag(U_temp)*ones(1, N_tx_sym);
U = U_temp./norm_mat;

Feedback_mat = U - eye(N_tx_sym);
% Feedback_mat = U;
Feedforward_mat = U*MMSE_mat;

data_ff = Feedforward_mat*data_rcv.';
decoded_soft = zeros(N_tx_sym, 1);
decoded_hard = zeros(N_tx_sym, 1);
decoded_soft(N_tx_sym) = data_ff(N_tx_sym);
decoded_hard(N_tx_sym) = modulator(demodulator(decoded_soft(N_tx_sym), S), S);
% decoded_hard(N_tx_sym) = modulator(demodulator(decoded_soft(N_tx_sym), mod_level), mod_level);

for k = N_tx_sym-1:-1:1
   decoded_soft(k) = data_ff(k) - Feedback_mat(k, :)*decoded_hard;
   data_demod = demodulator(decoded_soft(k), S);
   decoded_hard(k) = modulator(data_demod, S);
   % data_demod = demodulator(decoded_soft(k), mod_level);
   % decoded_hard(k) = modulator(data_demod, mod_level);
end


