function [x_hat,iters,t_RXiter,t_RXfull] = equalizer_CMC_MMSE_AWGN(y_tilde,H_tilde,N,M,Lp,Ln,Es,N0,S,N_iters,known_mask,known_x)
% BUG FIX (2026-09-14): added known_mask/known_x (optional, backward-
% compatible - omitting them reproduces the exact pre-fix behavior).
% Found via a direct paired comparison against equalizer_CMC_MMSE_native.m
% (a fresh, independent implementation of the same CP-Free ODDM.pdf
% Eq. 19-27 algorithm): this function had no mechanism to mark any delay
% layer n as KNOWN (e.g. a zero-padding guard) - EVERY layer, guard
% included, was hard-sliced onto the QPSK alphabet S every iteration
% (line ~78-80 below) and that wrong, nonzero slice was then fed back as
% ISI for neighboring layers within L1/L2 reach of the guard boundary.
% Same root cause, same fix pattern, as the equalizer_SIC_MMSE.m bug
% documented elsewhere in this project (hard-slicing a known-zero layer
% corrupts SIC feedback for nearby data layers) - confirmed empirically:
% a 300-frame paired BER test (v=500, EbN0=16dB) found this function gave
% 2.37x the BER of the known-mask-respecting equalizer_CMC_MMSE_native.m,
% despite the two being otherwise mathematically identical (verified:
% data-layer x_hat matched EXACTLY, to 0.0000e+00, before this fix, on an
% independent smaller sample - the corruption is real but probabilistic,
% only showing up in aggregate over enough frames).

% Start runtime
tStartRX = tic;

if nargin < 12 || isempty(known_mask)
    known_mask = false(M,1);
    known_x = zeros(N*M,1);
end
known_x_masked = known_x;
known_x_masked(~(repelem(known_mask(:),N) > 0)) = 0;

% Find all possible Lambda_n matrices and Theta_n matrices, and the MMSE
% matrix W_n they produce. W_n depends only on Lambda_n, N0 and Es - none
% of which change across equalizer iterations - so it is computed once
% here per block n and reused every iteration, instead of being
% recomputed (with a redundant pinv call) on every pass below.
possible_W_n = zeros(N*M,N);
for n = 0:M-1
    Lambda_n = zeros(N);
    for l = Ln:Lp
        selected_H_block1 = H_tilde((mod(n+l,M)*N)+1:(mod(n+l,M)+1)*N,(n*N)+1:(n+1)*N);
        Lambda_n_add = selected_H_block1' * selected_H_block1;
        Lambda_n = Lambda_n + Lambda_n_add;
    end
    possible_W_n((n*N)+1:(n+1)*N,:) = Lambda_n' * pinv(Lambda_n*Lambda_n' + (N0/Es)*Lambda_n);
end

% Predefine variables and start iterator equalizer. Known (e.g.
% zero-padding guard) layers are seeded with their TRUE value here and
% are NEVER re-touched below (skipped via known_mask(n+1) in the n-loop),
% so every ISI read of x_hat for a known interferer m automatically sees
% its correct, fixed value throughout every iteration - no separate
% known-value lookup needed in the ISI loop itself.
x_hat = known_x_masked;
flag_detector = true;
iters = 0;
iter_runtimes = [];
while iters < N_iters && flag_detector
    iters = iters + 1;

    % Start runtime
    tStartRXiter = tic;

    % Sweep through all M blocks of y_tilde (size Nx1)
    for n = 0:M-1
        if known_mask(n+1)
            continue;   % known layer - never re-detected, x_hat already holds its true value
        end
        % Create gamma_n with both for loops
        gamma_n = zeros(N,1);

        for l = Ln:Lp
            % Calculate ISI from chosen x_hat_l
            ISI = zeros(N,1);
            for k = Ln:Lp
                if k ~= l
                    % Update block indices
                    index1 = mod(n+l,M);
                    index2 = mod(n+l-k,M);

                    % Select current blocks and make ISI to add
                    selected_H_block1 = H_tilde((index1*N)+1:(index1+1)*N,(index2*N)+1:(index2+1)*N);
                    selected_x_block = x_hat((index2*N)+1:(index2+1)*N);
                    ISI_add = selected_H_block1 * selected_x_block;

                    % Add to ISI
                    ISI = ISI + ISI_add;
                end
            end
            % Update block indices
            index1 = mod(n+l,M);
            index2 = n;

            % Create y_tilde for current l
            y_tilde_l = y_tilde((index1*N)+1:(index1+1)*N);

            % Create y_hat
            y_hat_l = y_tilde_l - ISI;

            % Add to Gamma_n
            selected_H_block1 = H_tilde((index1*N)+1:(index1+1)*N,(index2*N)+1:(index2+1)*N);
            gamma_n = gamma_n + selected_H_block1' * y_hat_l;
        end

        % Select the precomputed MMSE matrix for this block (iteration-invariant)
        W_n = possible_W_n((n*N)+1:(n+1)*N,:);

        % Create x_hat for current block and push to stack
        x_hat_n = W_n * gamma_n;

        % Hard detection for block
        dist1 = abs(x_hat_n.' - S).^2;
        [~,min_index1] = min(dist1);
        x_hat_n = S(min_index1);

        % Map hard encoded x_hat_n to x_hat
        x_hat((index2*N)+1:(index2+1)*N) = x_hat_n;
    end

    % Check if duplicate result is found, break if true
    if iters > 1
        if last_x_hat == x_hat
            flag_detector = false;
        end
    end

    % Save last x_hat
    last_x_hat = x_hat;

    % Stop runtime
    t_RXiter = toc(tStartRXiter);
    iter_runtimes = [iter_runtimes t_RXiter]; %#ok<AGROW>

end

% Stop runtime
t_RXiter = mean(iter_runtimes);
t_RXfull = toc(tStartRX);