function [x_hat,iter,t_RXiter,t_RXfull] = equalizer_SIC_MMSE(r,G,N,M,L,Es,N0,S_alphabet,N_iters)
% Time-domain SIC-MMSE receiver seen in:
% "Iterative MMSE Detection for Orthogonal Time Frequency Space Modulation"
%     by Dr. Jinhong Yuan and Dr. Hai Lin
%
% Coded by JRW, 1/21/2026
% Ported verbatim from MUSIC OTFS Channel Estimation project, 9/2026

% Start runtime
tStartRX = tic;

% Initialize variables
tol = 1e-5;
s_hat = zeros(M,N,N_iters);
r_block = reshape(r,M,N);
F_N = gen_DFT(N);

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

            % Select sub-channel for (k+1)-th layer, and its first column
            H_e = G_n((1+k):(L+1+k),(1+k):(L+1+k));
            g = H_e(:,1);

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
            end

            % Generate MMSE matrix for i-th element
            w_MMSE = g' * pinv(H_e * H_e' + N0/Es * eye(L+1));

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
