function decoded_hard_shuffle = equalizer_MP()

%************************************
% block message passing (BMP)
% threshold for discarding elements in H_mat, in dB
threshold = 1e-6;
variable_idx_vec = [1:NM];
H_mat_logical = (abs(H_mat_DD_shuffle) > threshold);

for iM = 0:M-1
    H_mat_logical_sub_wide = H_mat_logical(iM*N+[1:N], :);
    % with block MP, some of the columns of H_mat_logical might be all zero. We need to mask out those columns
    variable_node_mask_mat(iM+1, :) = (sum(H_mat_logical_sub_wide, 1)>0);
end


N_sub_iteration = 8;
% initialization of data_tx_soft
data_tx_soft = 1/mod_level*ones(mod_level, NM);
% used to determine whether to update the soft information of a given variable
% A variable will only be updated if the convergence_vec > pre_convergence_vec
convergence_vec = max(data_tx_soft, [], 1);


N_block_itr = 1;
for i_block_itr = 1:N_block_itr
    for iM = 0:M-1
        data_vec_sub = data_rx_DD_vec_shuffle(iM*N+[1:N]);
        current_variable_node_mask_vec = variable_node_mask_mat(iM+1, :);
        prev_convergence_vec = convergence_vec(current_variable_node_mask_vec);

        H_mat_sub_wide = H_mat_DD_shuffle(iM*N+[1:N], current_variable_node_mask_vec);


        current_data_tx_soft = data_tx_soft(:, current_variable_node_mask_vec);


        [decoded_hard_sub, decoded_soft_sub] = message_passing_soft(data_vec_sub, H_mat_sub_wide, current_data_tx_soft, mod_sym_vec, N_sub_iteration, noise_var);

        % update data_tx_soft
        current_convergence_vec = max(decoded_soft_sub, [], 1);
        update_flag = (current_convergence_vec > prev_convergence_vec);
        current_idx_vec = variable_idx_vec(current_variable_node_mask_vec);
        update_idx_vec = current_idx_vec(update_flag);
        data_tx_soft(:, update_idx_vec) = decoded_soft_sub(:, update_flag);

        convergence_vec(update_idx_vec) = current_convergence_vec(update_flag);

        if and(iM == 0, mm == 37)
            data_mod_shuffle(1:2)
            H_mat_logical
            decoded_soft_sub
            current_variable_node_mask_vec
            convergence_vec
        end

    end
end

% final decision
[max_val, max_idx] = max(data_tx_soft, [], 1);
decoded_hard_shuffle = mod_sym_vec(max_idx);
%******************************************