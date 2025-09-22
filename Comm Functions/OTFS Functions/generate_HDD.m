function H_DD = generate_HDD(obj,ambig_vals,ambig_t_range,ambig_f_range,L)

% Import settings from object
Fc = obj.Fc;
v_vel = obj.v_vel;
sbcar_spacing = obj.sbcar_spacing;
% select_filter = obj.select_filter;
N_tsyms = obj.N_tsyms;
M_sbcars = obj.M_sbcars;
T = 1 / sbcar_spacing;

% Create channel gain, delay and doppler shift values for each path
[chn_g,chn_tau,chn_v] = channel_generation(Fc,v_vel);

% % Create each element of H_TF
% H_TF = zeros(N_tsyms*M_sbcars);
% for n1 = 0:N_tsyms-1
%     for n2 = 0:N_tsyms-1
%         % Compute time delay term
%         time_diff = (n1 - n2) * T;
%         time_diff_terms = time_diff - chn_tau;
% 
%         % Compute second exponential values
%         exp_vals2 = exp(1j*2*pi*chn_v*n2*T);
% 
%         for m1 = 0:M_sbcars-1
%             for m2 = 0:M_sbcars-1
%                 % Compute frequency delay term
%                 sbcar_diff = (m1 - m2) * sbcar_spacing;
%                 sbcar_diff_terms = sbcar_diff - chn_v;
% 
%                 % Decide if element is within render scope
%                 if abs(m1-m2) <= L
% 
%                     % Compute ambiguity function values
%                     xambig_vals = interp2( ...
%                         ambig_f_range, ...   % Doppler/freq axis
%                         ambig_t_range, ...   % Delay/time axis
%                         ambig_vals, ...      % Ambiguity surface
%                         sbcar_diff_terms, ...% Query freq (vector over paths)
%                         time_diff_terms, ... % Query delay (vector over paths)
%                         'linear', 0);        % Outside grid = 0
% 
%                     % Compute first exponential values
%                     exp_vals1 = exp(1j.*2.*pi.*(chn_v + m2.*sbcar_spacing).*(time_diff_terms));
% 
%                     % Calculate h value
%                     h = sum(chn_g .* xambig_vals .* exp_vals1 .* exp_vals2);
% 
%                     % Set h to be a value in H_TF
%                     H_TF(m1*N_tsyms+n1+1,m2*N_tsyms+n2+1) = h;
% 
%                 end
%             end
%         end
%     end
% end

% --- Vectorized block-wise construction of H_TF ---
N = N_tsyms;
M = M_sbcars;

% Precompute n-grid and flattened quantities used for every block
[n1_mat, n2_mat] = ndgrid(0:N-1, 0:N-1);   % n1 rows, n2 columns
time_diff = (n1_mat - n2_mat) * T;        % N x N
time_flat = time_diff(:);                 % (N^2) x 1
n2_flat = n2_mat(:);                      % (N^2) x 1, needed for exp2

% Row vectors of channel path params (1 x P) for broadcasting
chn_g_row   = chn_g(:).';   % 1 x P
chn_v_row   = chn_v(:).';   % 1 x P
chn_tau_row = chn_tau(:).'; % 1 x P

% Prepare H_TF
H_TF = zeros(N*M, N*M);

% Loop over blocks (m1,m2) but vectorize all n1,n2/path computations inside
for m1 = 0:M-1
    rowStart = m1 * N + 1;
    rowEnd   = (m1+1) * N;
    for m2 = 0:M-1
        colStart = m2 * N + 1;
        colEnd   = (m2+1) * N;

        if abs(m1 - m2) <= L
            % scalar sbcar diff for this block
            sbcar_diff = (m1 - m2) * sbcar_spacing;

            % Build query arrays for interp2:
            % sbcar_query: (N^2) x P  (frequency queries for each path)
            % time_query:  (N^2) x P  (delay queries for each path)
            sbcar_query = repmat(sbcar_diff - chn_v_row, N^2, 1);              % (N^2 x P)
            time_query  = bsxfun(@minus, time_flat, chn_tau_row);            % (N^2 x P)

            % Interpolate ambiguity table for all query points and all paths
            % ambig_vals: rows -> ambig_t_range (time), cols -> ambig_f_range (freq)
            xambig_flat = interp2( ...
                ambig_f_range, ...   % x axis: frequency
                ambig_t_range, ...   % y axis: time
                ambig_vals, ...      % surface
                sbcar_query, ...     % query X (N^2 x P)
                time_query, ...      % query Y (N^2 x P)
                'linear', 0);        % outside grid -> 0
            % xambig_flat is (N^2) x P

            % Compute exponentials (all (N^2 x P) arrays)
            % exp_vals1 = exp(1j*2*pi*(chn_v + m2*sbcar_spacing) .* (time_diff_terms))
            exp_vals1_flat = exp(1j * 2*pi * bsxfun(@times, (chn_v_row + m2*sbcar_spacing), time_query));

            % exp_vals2 = exp(1j*2*pi*chn_v * n2 * T)  (depends only on n2 and chn_v)
            exp_vals2_flat = exp(1j * 2*pi * bsxfun(@times, (n2_flat * T), chn_v_row));  % (N^2 x P)

            % Multiply everything and sum over paths (dimension 2)
            product = bsxfun(@times, chn_g_row, xambig_flat .* exp_vals1_flat .* exp_vals2_flat); % (N^2 x P)
            h_flat = sum(product, 2);   % (N^2 x 1)

            % Reshape back to N x N block (columns correspond to n2, rows to n1)
            H_TF(rowStart:rowEnd, colStart:colEnd) = reshape(h_flat, N, N);

        else
            % outside render window -> leave block as zeros (already zero)
        end
    end
end
% --- end vectorized block construction ---

% Create DFT matrices and needed Kronecker product
F_N = generate_DFT(N_tsyms);
F_M = generate_DFT(M_sbcars);
F_cur = kron(F_M,F_N');

% Create final result
H_DD = F_cur' * H_TF * F_cur;