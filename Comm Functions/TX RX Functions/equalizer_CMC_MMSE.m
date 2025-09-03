function [x_hat,iters,t_RXiter,t_RXfull] = equalizer_CMC_MMSE(y_tilde,H_tilde,N,M,Lp,Ln,Es,N0,S,N_iters,R)

% Start runtime
tStartRX = tic;

% Find all possible Lambda_n matrices and Theta_n matrices
possible_Lambda_n = zeros(N*M,N);
possible_Theta_n = zeros(N*M,N);
for n = 0:M-1
    Lambda_n = zeros(N);
    Theta_n = zeros(N);
    for l = Ln:Lp
        selected_H_block1 = H_tilde((mod(n+l,M)*N)+1:(mod(n+l,M)+1)*N,(n*N)+1:(n+1)*N);
        Lambda_n_add = selected_H_block1' * selected_H_block1;
        Lambda_n = Lambda_n + Lambda_n_add;

        for g = Ln:Lp
            selected_H_block2 = H_tilde((mod(n+g,M)*N)+1:(mod(n+g,M)+1)*N,(n*N)+1:(n+1)*N);
            selected_R_block = R((mod(n+l,M)*N)+1:(mod(n+l,M)+1)*N,(mod(n+g,M)*N)+1:(mod(n+g,M)+1)*N);
            Theta_n_add = selected_H_block1' * selected_R_block * selected_H_block2;
            Theta_n = Theta_n + Theta_n_add;
        end

    end
    possible_Lambda_n((n*N)+1:(n+1)*N,:) = Lambda_n;
    possible_Theta_n((n*N)+1:(n+1)*N,:) = Theta_n;
end

% Predefine variables and start iterator equalizer
x_hat = zeros(N*M,1);
flag_detector = true;
iters = 0;
iter_runtimes = [];
while iters < N_iters && flag_detector
    iters = iters + 1;

    % Start runtime
    tStartRXiter = tic;

    % Sweep through all M blocks of y_tilde (size Nx1)
    for n = 0:M-1
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

        % Select Lambda_n from possible matrices
        Lambda_n = possible_Lambda_n((n*N)+1:(n+1)*N,:);
        Theta_n = possible_Theta_n((n*N)+1:(n+1)*N,:);

        % Create MMSE matrix for iterative solver
        W_n = Lambda_n' / (Lambda_n*Lambda_n' + (N0/Es)*Theta_n);

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