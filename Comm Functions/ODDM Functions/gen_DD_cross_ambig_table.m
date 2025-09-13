function [ambig_vals,ambig_t_range,ambig_f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,v,shape,alpha,Q,res)

% Add redundancy for rectangular pulses
if shape == "rect"
    Q = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Define file and folder names
file_name = sprintf('DD_table_N%d_M%d_T%d,Fc%d_v%d_%s_a%.2f_Q%d_%dres',N,M,T,Fc,v,shape,alpha,Q,res);
folder_name = 'Pre-rendered Lookup Tables/ODDM DD Cross-Ambiguity Tables';
full_path = folder_name + "/" + file_name + ".mat";

% Check if the folder exists, and create it if it does not
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% Define ambiguity table
try
    
    % Load file if it exists
    load(full_path)

catch

    % fprintf("File not found for delay-Doppler cross ambiguity table...\nGenerating...\n")

    % Parameters
    Ts = T / M;
    F0 = 1 / (N*T);

    % Generate channel for max tau
    [~,tau_i,~] = channel_generation(Fc,v);
    max_tau = max(tau_i);

    % Set ambiguity range
    dt = Ts / res;
    df = F0 / res;
    ambig_t_lim_lower = M*Ts + max_tau;
    ambig_t_lim_upper = M*Ts;
    ambig_f_lim = (M-1)*F0 + ((v * (1000/3600))*Fc) / (physconst('LightSpeed'));
    ambig_t_lim_lower = ceil(ambig_t_lim_lower / dt) * dt;
    ambig_t_lim_upper = ceil(ambig_t_lim_upper / dt) * dt;
    ambig_f_lim = ceil(ambig_f_lim / df) * df;
    ambig_t_range = -ambig_t_lim_lower:dt:ambig_t_lim_upper;
    ambig_f_range = -ambig_f_lim:df:ambig_f_lim;

    % Render all ambiguity values
    ambig_vals = zeros(res);
    % update_vals = floor(length(ambig_t_range)*(linspace(.01,.99,34)));
    % fprintf("|           Progress:            |\n")
    for k = 1:length(ambig_t_range)
        % if ismember(k,update_vals)
        %     fprintf("x")
        % end
        for l = 1:length(ambig_f_range)
            ambig_vals(k,l) = DD_cross_ambig(ambig_t_range(k),ambig_f_range(l),N,M,T,shape,alpha,Q,res);
        end
    end
    % fprintf("\nComplete!\n\n")

    % Save to be used later
    save(full_path,"ambig_vals","ambig_t_range","ambig_f_range")
end

