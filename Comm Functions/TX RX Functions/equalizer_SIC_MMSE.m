function [x_hat,iter,t_RXiter,t_RXfull] = equalizer_SIC_MMSE(r,G,N,M,L,Es,N0,S_alphabet,N_iters,G_full)
% Time-domain SIC-MMSE receiver seen in:
% "Iterative MMSE Detection for Orthogonal Time Frequency Space Modulation"
%     by Dr. Jinhong Yuan and Dr. Hai Lin
%
% Coded by JRW, 1/21/2026
% Ported verbatim from MUSIC OTFS Channel Estimation project, 9/2026
% Extended 9/2026 (CWS-side only) to also cancel cross-time-symbol-block
% leakage using G_full: unlike OTFS, ODDM/CP-Free ODDM's channel is not
% exactly block-diagonal per time symbol, so a small amount of energy
% genuinely bleeds into the immediately adjacent block on either side
% (verified >99% of all off-block leakage energy is nearest-neighbor, both
% in this project and the sibling ODDM estimation paper project - see
% notes/AMBIGUITY_TABLE_AUDIT.md, "2026-09-09 follow-up"). Cancelled using each
% neighbor's own previous-iteration full-block estimate (zero on iter=1) -
% empirically equivalent to a same-iteration-where-available "mixed"
% timing scheme, per that same investigation, so the simpler scheme is
% used here.

% Start runtime
tStartRX = tic;

% Initialize variables
tol = 1e-5;
s_hat = zeros(M,N,N_iters);
r_block = reshape(r,M,N);
F_N = gen_DFT(N);

% Precompute the iteration-invariant MMSE row vector w_MMSE for every
% (k,n) layer/time-symbol pair. w_MMSE depends only on the channel (via
% G_n/H_e/g) and N0/Es - never on the SIC iteration or the ISI-cancelled
% residual - so it is computed once here instead of being recomputed
% (with a redundant pinv call) on every iteration below.
possible_w_MMSE = zeros(M-L,N,L+1);
for k = 0:M-L-1
    for n = 0:N-1
        G_n = G((n*M+1):((n+1)*M),(n*M+1):((n+1)*M));
        H_e = G_n((1+k):(L+1+k),(1+k):(L+1+k));
        g = H_e(:,1);
        possible_w_MMSE(k+1,n+1,:) = g' * pinv(H_e * H_e' + N0/Es * eye(L+1));
    end
end

% Precompute the cross-block leakage matrices: rows k..k+L of block n
% against the FULL M columns of its immediate neighbor blocks n-1,n+1
% (circular in the time-symbol/Doppler index, matching the periodic
% Kronecker-DFT construction in sim_fun_ODDM_SIC_MMSE.m). Like w_MMSE,
% these depend only on the channel, never on iteration/data, so they're
% computed once here too. G_full is optional (test/backward-compat use
% only - the production call site always supplies it) so this can be
% skipped entirely when omitted.
do_cross_block = nargin >= 10 && ~isempty(G_full);
if do_cross_block
    G_prev = cell(M-L,N);
    G_next = cell(M-L,N);
    for k = 0:M-L-1
        for n = 0:N-1
            rows_global = (n*M+k+1):(n*M+L+1+k);
            n_prev = mod(n-1,N);
            n_next = mod(n+1,N);
            G_prev{k+1,n+1} = G_full(rows_global,(n_prev*M+1):((n_prev+1)*M));
            G_next{k+1,n+1} = G_full(rows_global,(n_next*M+1):((n_next+1)*M));
        end
    end
end

% Loop through each iteration
iter_runtimes = [];
for iter = 1:N_iters

    % Start runtime
    tStartRXiter = tic;

    % Loop through each layer
    for k = 0:M-L-1

        % Loop through each time symbol
        for n = 0:N-1

            % Select current layers block to equalize
            G_n = G((n*M+1):((n+1)*M),(n*M+1):((n+1)*M));

            % Select L+1 received elements
            r_n = r_block((k+1):(L+1+k),n+1);

            % Perform interference cancelation - complete after first loop
            r_n_tilde = r_n;
            for m = 0:k-1 % Remove ISI using this iteration's estimate
                g_m = G_n((k+1):(k+L+1),m+1);
                r_n_tilde = r_n_tilde - g_m * s_hat(m+1,n+1,iter);
            end
            if iter > 1 % Remove ISI using last iteration's estimate
                for m = k+1:k+L
                    g_m = G_n((k+1):(k+L+1),m+1);
                    r_n_tilde = r_n_tilde - g_m * s_hat(m+1,n+1,iter-1);
                end

                % Remove cross-block leakage from the immediate neighbor
                % time-symbols using their own previous-iteration estimate
                if do_cross_block
                    n_prev = mod(n-1,N);
                    n_next = mod(n+1,N);
                    r_n_tilde = r_n_tilde - G_prev{k+1,n+1} * s_hat(:,n_prev+1,iter-1);
                    r_n_tilde = r_n_tilde - G_next{k+1,n+1} * s_hat(:,n_next+1,iter-1);
                end
            end

            % Retrieve the precomputed MMSE row vector for this layer/time-symbol
            w_MMSE = reshape(possible_w_MMSE(k+1,n+1,:), 1, L+1);

            % Do soft equalization for s
            s_hat(k+1,n+1,iter) = w_MMSE * r_n_tilde;

        end

        % Get estimated DD symbols for all time symbols, one layer at a time
        x_DD_tildem = F_N * s_hat(k+1,:,iter).';

        % Perform hard detection and reassign s_hat
        costs = abs(S_alphabet.' - x_DD_tildem).^2;
        [~,idx] = min(costs,[],2);
        x_DD_hatm = S_alphabet(idx);
        s_hat(k+1,:,iter) = (F_N' * x_DD_hatm).';

    end

    % Check if should stop before N_iters is completed
    if iter > 1
        if norm(s_hat_last - s_hat(:,:,iter)) < tol
            break;
        end
    end
    s_hat_last = s_hat(:,:,iter);

    % Stop runtime
    t_RXiter = toc(tStartRXiter);
    iter_runtimes = [iter_runtimes t_RXiter]; %#ok<AGROW>

end

% Export results
X_hat = s_hat_last * F_N;
x_hat = X_hat(:);

% Stop runtime
t_RXiter = mean(iter_runtimes);
t_RXfull = toc(tStartRX);
