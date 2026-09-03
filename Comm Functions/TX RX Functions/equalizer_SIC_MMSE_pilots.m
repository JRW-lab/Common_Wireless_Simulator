function [x_hat,iter,t_RXiter,t_RXfull] = equalizer_SIC_MMSE_pilots(r,G,N,M,num_pilots,L,Es,N0,S_alphabet,N_iters)
% Time-domain SIC-MMSE receiver seen in:
% "Iterative MMSE Detection for Orthogonal Time Frequency Space Modulation"
%     by Dr. Jinhong Yuan and Dr. Hai Lin
%
% This code has been modified to assume that instead of N blocks we have P
% pilot symbols evenly spaced. Each pilot block is equalized independently
% and once a full time-symbol (row) estimate is found, then hard
% equalization happens. It is assumed that N*M/P is an integer.
%
% Coded by JRW, 2/5/2026

% Start runtime
tStartRX = tic;

% Pilot based parameters
pilot_pattern_size = 2*L+1;
block_size = N*M/num_pilots + 2*L + 1;
pilots_per_tsym = num_pilots / N;

% Initialize variables
tol = 1e-5;
s_hat = zeros(M,N,N_iters);
r_block = reshape(r,block_size,num_pilots);
F_N = gen_DFT(N);

% Set range of block taps
sbcar_vec = (0:block_size-pilot_pattern_size-1).';
if pilots_per_tsym < 1
    sbcar_mat = reshape(sbcar_vec,length(sbcar_vec)*pilots_per_tsym,1/pilots_per_tsym).';
    sbcar_vec = sbcar_mat(:);
end

% Precompute the iteration-invariant MMSE row vector w_MMSE for every
% (k,n) layer/time-symbol pair. w_MMSE depends only on the channel (via
% G_n/H_e/g) and N0/Es - never on the SIC iteration or the ISI-cancelled
% residual - so it is computed once here instead of being recomputed
% (with a redundant pinv call) on every iteration below.
num_k = block_size - pilot_pattern_size;
possible_w_MMSE = zeros(num_k,num_pilots,L+1);
for k = sbcar_vec.'
    for n = 0:num_pilots-1
        G_n = G((n*block_size+1):((n+1)*block_size),(n*block_size+1):((n+1)*block_size));
        H_e = G_n((1+k):(L+1+k),(1+k):(L+1+k));
        g = H_e(:,1);
        possible_w_MMSE(k+1,n+1,:) = g' * pinv(H_e * H_e' + N0/Es * eye(L+1));
    end
end

% Loop through each iteration
iter_runtimes = [];
for iter = 1:N_iters

    % Start runtime
    tStartRXiter = tic;

    % Loop through each layer
    for k = sbcar_vec.'

        % Loop through each time symbol
        for n = 0:num_pilots-1

            % Select current time symbol and subcarrier (true DD grid)
            if pilots_per_tsym >= 1
                tsym_sel = ceil((n+1)/pilots_per_tsym);
                sbcar_sel = k+1 + M*(mod(n,pilots_per_tsym))/pilots_per_tsym;
            else
                tsym_sel = 1 + n/pilots_per_tsym + (k > (M-1));
                sbcar_sel = mod(k,M)+1;
            end

            % Select current layers block to equalize
            G_n = G((n*block_size+1):((n+1)*block_size),(n*block_size+1):((n+1)*block_size));

            % Select L+1 received elements
            r_n = r_block((k+1):(L+1+k),n+1);

            % Perform interference cancelation - complete after first loop
            r_n_tilde = r_n;
            for m = 0:k-1
                % Translate to correct subcarrier
                if pilots_per_tsym >= 1
                    sbcar_sel2 = m+1 + M*(mod(n,pilots_per_tsym))/pilots_per_tsym;
                else
                    sbcar_sel2 = mod(m,M)+1;
                end

                % Remove ISI using this iteration's estimate
                g_m = G_n((k+1):(k+L+1),m+1);
                r_n_tilde = r_n_tilde - g_m * s_hat(sbcar_sel2,tsym_sel,iter);
            end
            if iter > 1
                for m = k+1:min(k+L,block_size-pilot_pattern_size-1)
                    % Translate to correct subcarrier
                    if pilots_per_tsym >= 1
                        sbcar_sel2 = m+1 + M*(mod(n,pilots_per_tsym))/pilots_per_tsym;
                    else
                        sbcar_sel2 = mod(k,M)+1;
                    end

                    % Remove ISI using last iteration's estimate
                    g_m = G_n((k+1):(k+L+1),m+1);
                    r_n_tilde = r_n_tilde - g_m * s_hat(sbcar_sel2,tsym_sel,iter);
                end
            end

            % Retrieve the precomputed MMSE row vector for this layer/time-symbol
            w_MMSE = reshape(possible_w_MMSE(k+1,n+1,:), 1, L+1);

            % Do soft equalization for s
            s_hat(sbcar_sel,tsym_sel,iter) = w_MMSE * r_n_tilde;

        end

        if pilots_per_tsym >= 1
            for i = 1:pilots_per_tsym
                % Get estimated DD symbols for all time symbols, one layer at a time
                sbcar_idx = k+1 + (i-1)*M/pilots_per_tsym;
                x_DD_tildem = F_N * s_hat(sbcar_idx,:,iter).';

                % Perform hard detection and reassign s_hat
                costs = abs(S_alphabet.' - x_DD_tildem).^2;
                [~,idx] = min(costs,[],2);
                x_DD_hatm = S_alphabet(idx);
                s_hat(sbcar_idx,:,iter) = (F_N' * x_DD_hatm).';
            end
        elseif k > (M-1)
            % Get estimated DD symbols for all time symbols, one layer at a time
            sbcar_idx = mod(k,M) + 1;
            x_DD_tildem = F_N * s_hat(sbcar_idx,:,iter).';

            % Perform hard detection and reassign s_hat
            costs = abs(S_alphabet.' - x_DD_tildem).^2;
            [~,idx] = min(costs,[],2);
            x_DD_hatm = S_alphabet(idx);
            s_hat(sbcar_idx,:,iter) = (F_N' * x_DD_hatm).';
        end

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
