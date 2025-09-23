function [decoded_hard_shuffle,iters,t_RXiter,t_RXfull] = equalizer_MP(data_rx_DD_vec_shuffle, H_mat_DD_shuffle, mod_level, M, N, mod_sym_vec, noise_var)

% Start runtime
tStartRXiter = tic;

% Needed variables
NM = M*N;
N_block_itr = 1;
N_sub_iteration = 8;

% WUJ
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

% initialization of data_tx_soft
data_tx_soft = 1/mod_level*ones(mod_level, NM);
% used to determine whether to update the soft information of a given variable
% A variable will only be updated if the convergence_vec > pre_convergence_vec
convergence_vec = max(data_tx_soft, [], 1);


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

        % if and(iM == 0, mm == 37)
        %     data_mod_shuffle(1:2)
        %     H_mat_logical
        %     decoded_soft_sub
        %     current_variable_node_mask_vec
        %     convergence_vec
        % end

    end
end

% final decision
[max_val, max_idx] = max(data_tx_soft, [], 1);
decoded_hard_shuffle = mod_sym_vec(max_idx);
%******************************************

% Stop runtime
t_RXiter = toc(tStartRXiter);
t_RXfull = t_RXiter;
iters = N_sub_iteration;

end



%% SUBFUNCTIONS
function [data_detect, data_pmf] = message_passing_soft(data_rx, H_mat, data_tx_soft, mod_symbols, N_iteration, noise_var, damping_factor)
% function[data_detect, data_pmf] = messag_passing(data_rx, H_mat, data_tx_soft, mod_symbols, N_iteration, noise_var, damping_factor)
%
%	
% The code is based on 
% [1] P. Raviteja, K. T. Phan, Q. Jin, Y. Hong and E. Viterbo, "Low-complexity iterative detection for orthogonal time frequency space modulation," 2018 IEEE Wireless Communications and Networking Conference (WCNC), 2018, pp. 1-6, doi: 10.1109/WCNC.2018.8377159. 
%


%damping_factor = 1;
if nargin == 6
	damping_factor = 0.7;
end



% the convergence rate is used to terminate the iterations
convergence_rate_prev = -0.1;
convergene_rate = 0;

[N_check, N_variable] = size(H_mat);

if (length(data_rx) ~= N_check)
    error('data_rx and H_mat size mismatch.')
end

% data_rx needs to be a column vector
if isrow(data_rx)
    data_rx = data_rx.';
end

% mod_symbols needs to be a column vector
if isrow(mod_symbols)
    mod_symbols = mod_symbols.';
end



mod_level = length(mod_symbols);
abs_square_mod_symobls = abs(mod_symbols).^2;

threshold = 1e-6; 
H_mat_logical = (abs(H_mat) > threshold);

abs_square_H_mat = abs(H_mat).^2;

% the node degree (the number of variable nodes connnected to each variable node) for each check node
check_node_degree = sum(H_mat_logical, 2);
max_check_node_degree = max(check_node_degree);

% the node degree (the number of check nodes connnected to each variable node) for each variable node
variable_node_degree = sum(H_mat_logical, 1).';
max_variable_node_degree = max(variable_node_degree);

% The index matrix of all variable nodes connected to a check node
% Dimension: N_check x max_check_node_degree
% the k-th row contains the indices of all variable nodes connected to
% the k-th check node
check_node_idx_mat = zeros(N_check, max_check_node_degree);
check_node_chnn_mat = zeros(N_check, max_check_node_degree);
variable_idx_vec = [1:N_variable];
for iCheck = 1:N_check
    check_node_idx_mat(iCheck, [1:check_node_degree(iCheck)]) = variable_idx_vec(H_mat_logical(iCheck, :));
    check_node_chnn_mat(iCheck, [1:check_node_degree(iCheck)]) = H_mat(iCheck, H_mat_logical(iCheck, :));
    check_node_abs_square_chnn_mat(iCheck, [1:check_node_degree(iCheck)]) = abs_square_H_mat(iCheck, H_mat_logical(iCheck, :));
end



% The index matrix of all check nodes connected to a variable node
% Dimension: N_variable x max_variable_node_degree
% the k-th row contains the indices of all check nodes connected to
% the k-th variable node
variable_node_idx_mat = zeros(N_variable, max_variable_node_degree);
check_idx_vec = [1:N_check];
for iVariable = 1:N_variable
    variable_node_idx_mat(iVariable, [1:variable_node_degree(iVariable)]) = check_idx_vec(H_mat_logical(:, iVariable).');
end



  
% Initilization
% messages from variable nodes to check nodes
% Dimension: mod_level x N_variable x N_check
% the (:, n, k)-th column is the message from the n-th variable node to the k-th check node
% for all modulation constellation symbols
% the message will be initialized by using the soft input data_tx_soft (dimension mod_level x N_variable)
%message_var_to_check = 1/mod_level*ones(mod_level, N_variable, N_check);
for iCheck = 1:N_check
	message_var_to_check(:, :, iCheck) = data_tx_soft;
end



% messages (including both mean and variance) from check nodes to variable nodes
message_mean_check_to_var = zeros(N_check, N_variable);
message_variance_check_to_var = zeros(N_check, N_variable);


joint_pmf_mat = zeros(N_variable, mod_level);



for iIteration = 1:N_iteration

	% used to check stopping condition
	max_pmf_diff = 0;

    %-----------------------------------------------------------------
    % Update the messages from the check nodes to the variable nodes
    %-----------------------------------------------------------------
    for iCheck = 1:N_check
        current_check_node_degree = check_node_degree(iCheck);
        connected_var_nodes_id = check_node_idx_mat(iCheck, [1:current_check_node_degree]);
        

       	% the messages from all connected variable nodes
        prev_message_mat = message_var_to_check(:, connected_var_nodes_id, iCheck);
        
       	current_chnn_vec = check_node_chnn_mat(iCheck, [1:current_check_node_degree]);
        current_abs_square_chnn_vec = check_node_abs_square_chnn_mat(iCheck, [1:current_check_node_degree]);
        
        % calculate the mean based on eqn. (19) from [1] and the variance based on eqn. (20) from [1]
        current_mean_vec = zeros(1, current_check_node_degree);
        current_var_vec = zeros(1, current_check_node_degree);
        for iConnectedVar = 1:current_check_node_degree
        	current_connected_var_id = connected_var_nodes_id(iConnectedVar);
        	current_chnn = current_chnn_vec(iConnectedVar);
        	current_abs_square_chnn = current_abs_square_chnn_vec(iConnectedVar);
        	temp_mean = 0;
        	temp_var = 0;
        	for iMod = 1:mod_level
        		temp_mean = temp_mean + prev_message_mat(iMod, iConnectedVar)*mod_symbols(iMod);
        		temp_var = temp_var + prev_message_mat(iMod, iConnectedVar)*abs_square_mod_symobls(iMod);
        	end
        	current_mean = temp_mean*current_chnn;
        	current_var = (temp_var - abs(temp_mean)^2)*current_abs_square_chnn;
        	current_mean_vec(iConnectedVar) = current_mean;
        	current_var_vec(iConnectedVar) = current_var;
        end
        current_total_mean = sum(current_mean_vec);
        current_total_var = sum(current_var_vec) + noise_var;
        
        for iConnectedVar = 1:current_check_node_degree
        	current_var_node = connected_var_nodes_id(iConnectedVar);
        	current_mean = current_total_mean - current_mean_vec(iConnectedVar);
        	current_var = current_total_var - current_var_vec(iConnectedVar);
        	
        	
        	message_mean_check_to_var(iCheck, current_var_node) = current_mean;
        	message_variance_check_to_var(iCheck, current_var_node) = current_var;
        	
            
            
            
        end
	end
        
            
    %-----------------------------------------------------------------
    % Update the messages from the variable nodes to the check nodes
    %-----------------------------------------------------------------
    
    % first save the message from the previous iteration
    prev_message_var_to_check = message_var_to_check;
    
    for iVariable = 1:N_variable
        current_variable_node_degree = variable_node_degree(iVariable);
        connected_check_nodes_id = variable_node_idx_mat(iVariable, [1:current_variable_node_degree]);
        
        seq_check_node_idx_vec = [1:current_variable_node_degree];
        
        
        
        mean_vec = message_mean_check_to_var(connected_check_nodes_id, iVariable);
        variance_vec = message_variance_check_to_var(connected_check_nodes_id, iVariable);
        
        
        % calculate the overall log-likelihood function
        log_likelihood_mat = zeros(current_variable_node_degree, mod_level);
        
        
        for iMod = 1:mod_level
        	for iCheck = 1:current_variable_node_degree
        		current_check_node_id = connected_check_nodes_id(iCheck);
            	log_likelihood_mat(iCheck, iMod) = -abs(data_rx(current_check_node_id)-mean_vec(iCheck)-H_mat(current_check_node_id, iVariable)*mod_symbols(iMod)).^2./variance_vec(iCheck);
            end
        end
        
        %for iMod = 1:mod_level
        %    log_likelihood_mat(:, iMod) = -abs(data_rx(connected_check_nodes_id)-mean_vec-H_mat(connected_check_nodes_id, iVariable)*mod_symbols(iMod)).^2./variance_vec;
        %end
        
        
       % calculate the joint likelihood function of the current variable node
 	   joint_log_likelihood = sum(log_likelihood_mat, 1); 
 	   % normalized the likelihood function into a probability mass function (PMF). The normalization is performed in the 
       % log domain to avoid overflow
       % For example, 
       %    e(x)/[e(x)+e(y)] is performed as 
       %	 x - log[e(x)+e(y)] = x - logsum(x, y)
       %	
       %
       % perform sum in the log domain to avoid overflow
       joint_log_likelihood_sum = logsum(joint_log_likelihood);
       % perform normalization in the log domain
       joint_log_pmf = joint_log_likelihood - joint_log_likelihood_sum;
       joint_pmf_mat(iVariable, :) = exp(joint_log_pmf);

 
        
        
        for iCheck = 1:current_variable_node_degree
        	current_check_node_id = connected_check_nodes_id(iCheck);
        
        	% new code
        	current_joint_log_likelihood = joint_log_likelihood - log_likelihood_mat(iCheck, :);
        	% normalized the likelihood function into a probability mass function (PMF). The normalization is performed in the 
            % log domain to avoid overflow
            % For example, 
            %    e(x)/[e(x)+e(y)] is performed as 
            %	 x - log[e(x)+e(y)] = x - logsum(x, y)
            %	
            %
            % perform sum in the log domain to avoid overflow
            joint_log_likelihood_sum = logsum(current_joint_log_likelihood);
            % perform normalization in the log domain
            joint_log_pmf = current_joint_log_likelihood - joint_log_likelihood_sum;
            joint_pmf = exp(joint_log_pmf);
           
            prev_pmf = prev_message_var_to_check(:, iVariable, current_check_node_id).';          
             
            new_pmf = damping_factor*joint_pmf+(1-damping_factor)*prev_pmf;
            
            % debug
            if sum(isnan(new_pmf))
            	new_pmf
            	error('overflow');
            end

            
            % assign new message 
            message_var_to_check(:, iVariable, current_check_node_id) = new_pmf.';
           
        end 
    end
    
    
    % check whether we need to terminate the iterations
    convergence_rate =  sum(max(joint_pmf_mat,[],2)>0.99)/(N_variable);
    if convergence_rate==1
        joint_pmf_mat_final = joint_pmf_mat;
        break;
    elseif convergence_rate > convergence_rate_prev
        convergence_rate_prev = convergence_rate;
        joint_pmf_mat_final = joint_pmf_mat;
    elseif (convergence_rate < convergence_rate_prev - 0.2) && convergence_rate_prev > 0.95
        break;
    end
        
end


[max_val, max_idx] = max(joint_pmf_mat_final, [], 2);
data_detect = mod_symbols(max_idx);

data_pmf = joint_pmf_mat_final.';

end

function y = logsum(x)
% function y = logsum(x)
% 
% x is a vector in the log-domain
% y = log(exp(x1)+exp(x2)+...+exp(xN))
%
% algorithm:
% log(exp(x1) + exp(x2)) = log(exp(x1))+log(1+exp(x2-x1)) = x1 + log[1+exp(x2-x1)]
% the results are calculated recursively
%

x = real(x);

if length(x) == 1
	y = x;
elseif length(x) > 2
	temp_y = logsum(x(1:2));
	y = logsum([temp_y, x(3:end)]);
elseif length(x) == 2
	y = max(x) + log(1+exp(min(x)-max(x)));
else
	error('something is wrong!');
end
 
end