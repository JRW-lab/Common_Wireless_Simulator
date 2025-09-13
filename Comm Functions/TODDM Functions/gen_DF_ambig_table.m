function [Ambig_Table_Whole,t_range,phi_range] = gen_DF_ambig_table(U,N,M,T,Fc,vel,shape,alpha,Q,res,simplify)

% Add redundancy for rectangular and sinc pulses
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Parameters
Ts = T / M;
Fs = 1 / Ts;
[~, tau_i] = channel_generation(Fc,vel);

% Set ambiguity range
t_taps = (-M*Ts:Ts:M*Ts) - tau_i.';
max_Doppler = (vel * (1000/3600))*Fc / physconst('LightSpeed');
count = 0;
f_taps = zeros(2*U-1,res);
for u = -(U-1):(U-1)
    count = count + 1;
    f_taps(count,:) = linspace(u*Fs-max_Doppler,u*Fs+max_Doppler,res);
end
t_range = unique(sort(t_taps(:)));
phi_range = unique(sort(f_taps(:)));

% Filenames
folder_name = 'Pre-rendered Lookup Tables/TODDM DF Cross-Ambiguity Tables';
subfolder = fullfile(folder_name, ['/' char(sprintf('U%d_N%d_M%d_T%d_Fc%d_v%d_%s_a%.2f_Q%d_%dres_simplify%d',U,N,M,T,Fc,vel,shape,alpha,Q,res,simplify))]);
file_name = sprintf('U%d_N%d_M%d_T%d_Fc%d_v%d_%s_a%.2f_Q%d_%dres_simplify%d',U,N,M,T,Fc,vel,shape,alpha,Q,res,simplify);
main_path = subfolder + "/" + file_name + ".mat";

% % Begin progress bar
% if simplify
%     fprintf("Constructing simplified tri-orthogonal delay-Doppler cross ambiguity table...\n")
% else
%     fprintf("Constructing direct tri-orthogonal delay-Doppler cross ambiguity table...\n")
% end
% update_vals = floor(length(t_range)*length(phi_range)*N*N*(linspace(0.01,1,34)));
% fprintf("|           Progress:            |\n")
% count = 0;

try

    % Load file if it exists
    load(main_path)
    Ambig_Table_Whole(1);

    % % Increment progress bar
    % for i = 1:34
    %     fprintf("x")
    % end

catch

    Ambig_Table_Whole = cell(N);
    for n = 0:N-1
        for k = 0:N-1

            % Define file and folder names
            subfile_name = sprintf('DD_table_n%d_k%d_U%d_N%d_M%d_T%d_Fc%d_v%d_%s_a%.2f_Q%d_%dres_simplify%d',n,k,U,N,M,T,Fc,vel,shape,alpha,Q,res,simplify);
            sub_path = subfolder + "/" + subfile_name + ".mat";

            % Check if the folder exists, and create it if it does not
            if ~exist(folder_name, 'dir')
                mkdir(folder_name);
            end
            if ~exist(subfolder, 'dir')
                mkdir(subfolder);
            end

            % Define ambiguity table
            try

                % Load file if it exists
                load(sub_path)

                % % Increment progress bar
                % count = count + length(t_range)*length(phi_range);
                % for i = 1:length(update_vals)
                %     flag1 = update_vals(i) <= (count);
                %     flag2 = update_vals(i) > (count - length(t_range)*length(phi_range));
                %     if flag1 && flag2
                %         fprintf("x")
                %     end
                % end

            catch

                % Render all ambiguity values
                Ambig_Table = zeros(length(t_range),length(phi_range));

                for i = 1:length(t_range)
                    for j = 1:length(phi_range)

                        % % Increment progress bar
                        % count = count + 1;
                        % if ismember(count,update_vals)
                        %     fprintf("x")
                        % end

                        Ambig_Table(i,j) = Ank(n,k,t_range(i),phi_range(j),N,M,T,shape,alpha,Q,res,simplify);
                    end
                end

                % Save to be used later
                save(sub_path,"Ambig_Table","t_range","phi_range")
            end

            % Add ambiguity table to set
            Ambig_Table_Whole{n+1,k+1} = Ambig_Table;

        end

    end

    % % Add ambiguity set
    % save(main_path,"Ambig_Table_Whole","t_range","phi_range","-v7.3")

end

% fprintf("\nComplete!\n\n")