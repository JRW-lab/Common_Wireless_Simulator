function H_DD = generate_HDD(obj)

% Method 1
% Debugging inputs
% clear; clc;
% obj = comms_obj_OTFS;

% Import settings from object
Fc = obj.Fc;
v_vel = obj.v_vel;
sbcar_spacing = obj.sbcar_spacing;
select_filter = obj.select_filter;
N_tsyms = obj.N_tsyms;
M_sbcars = obj.M_sbcars;

% Initialize and fill in block matrix
% render_range = zeros(N_tsyms, N_tsyms*M_sbcars);
% for i = 1:N_tsyms/2:N_tsyms*M_sbcars
%     render_range(i:i+N_tsyms-1, i:i+N_tsyms-1) = ones(N_tsyms);
% end
render_range = ones(N_tsyms*M_sbcars);

% Internal settings
T = 1 / sbcar_spacing;

% Create channel gain, delay and doppler shift values for each path
[chn_g,chn_tau,chn_v] = channel_generation(Fc,v_vel);

% Create each element of H_TF
H_TF = zeros(N_tsyms*M_sbcars);
for n1 = 0:N_tsyms-1
    for n2 = 0:N_tsyms-1
        % Compute time delay term
        time_diff = (n1 - n2) * T;
        time_diff_terms = time_diff - chn_tau;

        % Compute second exponential values
        exp_vals2 = exp(1j*2*pi*chn_v*n2*T);

        for m1 = 0:M_sbcars-1
            for m2 = 0:M_sbcars-1
                % Compute frequency delay term
                sbcar_diff = (m1 - m2) * sbcar_spacing;
                sbcar_diff_terms = sbcar_diff - chn_v;

                % Decide if element is within render scope
                if render_range(m1*N_tsyms+n1+1, m2*N_tsyms+n2+1)
                % if true

                    % Compute ambiguity function values
                    if select_filter ~= "ideal"
                        xambig_vals = xambig(time_diff_terms, sbcar_diff_terms, T, select_filter);
                    else
                        xambig_vals = xambig(time_diff * ones(1,length(chn_v)), sbcar_diff * ones(1,length(chn_v)), T, select_filter);
                    end

                    % Compute first exponential values
                    exp_vals1 = exp(1j.*2.*pi.*(chn_v + m2.*sbcar_spacing).*(time_diff_terms));

                    % Calculate h value
                    h = sum(chn_g .* xambig_vals .* exp_vals1 .* exp_vals2);

                    % % Create sum to find resulting 'h'
                    % h1 = 0;
                    % for k = 1:length(chn_g)
                    %     if select_filter ~= "ideal"
                    %         xambig_val = xambig((n1-n2)*T - chn_tau(k),(m1-m2)*sbcar_spacing - chn_v(k),T,select_filter);
                    %     else
                    %         xambig_val = xambig((n1-n2)*T,(m1-m2)*sbcar_spacing,T,select_filter);
                    %     end
                    %     exp_val1 = exp(1j*2*pi*(chn_v(k) + m2*sbcar_spacing)*((n1-n2)*T - chn_tau(k)));
                    %     exp_val2 = exp(-1j*2*pi*chn_v(k)*n2*T);
                    %     h1 = h1 + chn_g(k) * xambig_val * exp_val1 * exp_val2;
                    % end

                    % Set h to be a value in H_TF
                    H_TF(m1*N_tsyms+n1+1,m2*N_tsyms+n2+1) = h;
                    % H_TF1(m1*N_tsyms+n1+1,m2*N_tsyms+n2+1) = h1;

                end
            end
        end
    end
end

% Create DFT matrices and needed Kronecker product
F_N = generate_DFT(N_tsyms);
F_M = generate_DFT(M_sbcars);
F_cur = kron(F_M,F_N');

% Create final result
H_DD = F_cur' * H_TF * F_cur;

% % Method 2
% % Create all possible combinations of n1,n2,m1,m2
% counter = 0;
% H_combos = zeros(N_tsyms^2 * M_sbcars^2,4);
% for n1 = 0:N_tsyms-1
%     for n2 = 0:N_tsyms-1
%         for m1 = 0:M_sbcars-1
%             for m2 = 0:M_sbcars-1
%                 counter = counter + 1;
%                 H_combos(counter,1) = n1;
%                 H_combos(counter,2) = n2;
%                 H_combos(counter,3) = m1;
%                 H_combos(counter,4) = m2;
%             end
%         end
%     end
% end
% 
% % Sweep through and evaluate every value
% for n1 = 0:N_tsyms-1
%     for n2 = 0:N_tsyms-1
%         for m1 = 0:M_sbcars-1
%             for m2 = 0:M_sbcars-1
%                 counter = counter + 1;
%                 H_combos(counter,1) = n1;
%                 H_combos(counter,2) = n2;
%                 H_combos(counter,3) = m1;
%                 H_combos(counter,4) = m2;
%             end
%         end
%     end
% end