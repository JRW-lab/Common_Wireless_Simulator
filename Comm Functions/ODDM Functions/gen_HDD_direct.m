function [H,L1,L2,Phi_i,tau_i,v_i] = gen_HDD_direct(T,N,M,Fc,v,Q,Ambig_Table,t_offset,CP)
% This generates the delay-Doppler channel matrix for an ODDM system
% INPUTS:
%   T:              Symbol interval
%   N:              number of time symbols
%   M:              number of subcarriers
%   Fc:             carrier frequency
%   v:              maximum vehicle velocity
%   Q:              range of elementary pulse defined as [0,Q*Ts]
%   Ambig_Table:    prerendered cross-ambiguity table
%
% Coded by Jeremiah Rhys Wimer, 2/27/2025

% Import values from ambiguity table
ambig_vals = Ambig_Table.vals;
ambig_t_range = Ambig_Table.t_range;
ambig_f_range = Ambig_Table.f_range;

% Parameters
Ts = T / M;
F0 = 1 / (N*T);

% Generate channel
[Phi_i,tau_i,v_i] = channel_generation(Fc,v);

% Define channel taps
L1 = Q + 1;
L2 = Q + 1 + floor(max(tau_i) / Ts);

% Set up L tap range
range1 = -M + (1:L2);
range2 = -L1:L2;
range3 = M - (L1:-1:1);
L_range = [range1, range2, range3].';

% Find all combinations of m,n,l,k - then only select valid ones
if CP
    M_cp = L1+L2;
else
    M_cp = 0;
end
m_range = -M_cp:(M-1);
l_range = -M_cp:(M-1);
n_range = 0:(N-1);
k_range = 0:(N-1);
tap_combos = combvec(m_range,n_range,l_range,k_range).';
l_minus_m = tap_combos(:,3) - tap_combos(:,1);

needed_indices = ismember(l_minus_m,L_range);
needed_combos = tap_combos(needed_indices,:);
m = needed_combos(:,1);
n = needed_combos(:,2);
l = needed_combos(:,3);
k = needed_combos(:,4);

% Find indices in ambiguity values corresponding with the coefficient eq.
ambig_t_locs = (l-m).*Ts - tau_i + t_offset;
ambig_f_locs = (k-n).*F0 - v_i;
% ambig_t_indices = knnsearch(ambig_t_range(:), ambig_t_locs(:));
% ambig_t_indices = reshape(ambig_t_indices, size(ambig_t_locs));
% ambig_f_indices = knnsearch(ambig_f_range(:), ambig_f_locs(:));
% ambig_f_indices = reshape(ambig_f_indices, size(ambig_f_locs));
% 
% % Convert subscripts to linear indices and select ambiguity values
% linear_indices = sub2ind(size(ambig_vals), ambig_t_indices, ambig_f_indices);
% ambig_inst = ambig_vals(linear_indices);

%%% More efficient way to find ambiguity values
ambig_inst = interp2(ambig_f_range, ambig_t_range, ambig_vals, ...
                     ambig_f_locs, ambig_t_locs, 'linear', 0);

% Define all elements of h summation, then sum across 2nd dimension
h_sum = exp(-1j.*2.*pi.*n.*m./(N.*M)) .* Phi_i .* exp(1j.*2.*pi.*(v_i + n.*F0).*(l.*Ts-tau_i+t_offset)) .* ambig_inst;
h = sum(h_sum,2);

% Find all linear indices for H matrix
H = zeros(length(unique(m))*length(unique(n)));
H_index1 = l*N + k+1 - min(m_range)*N;
H_index2 = m*N + n+1 - min(m_range)*N;
linear_indices = sub2ind(size(H), H_index1, H_index2);

% Assign values from h to the corresponding locations in H
H(linear_indices) = h;

% % Delete unneeded blocks if using CP
% if CP
%     H(1:(size(H,1)-M*N-L1-L2),:) = 0;
% end

1;

% % ---- Map to Toeplitz matrix with CP ----
% % Total length including CP
% blk_len = M*N + L1 + L2;  
% H = zeros(blk_len, blk_len);
% 
% % Offset for CP in the output
% cp_offset = L1;
% 
% for idx = 1:length(h)
%     % Output index (with CP offset)
%     row = cp_offset + l(idx)*N + k(idx) + 1;
% 
%     % Input index
%     col = m(idx)*N + n(idx) + 1;
% 
%     % Ensure row/col within block length
%     if row <= blk_len && col <= blk_len
%         H(row, col) = H(row, col) + h(idx);
%     end
% end
