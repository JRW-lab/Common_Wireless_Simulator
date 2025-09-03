function [HDD,L1,L2] = gen_HDD_triorth(U,N,M,T,Q,Fc,vel,Ambig_Table,t_offset)

% Parameters
Ts = T / M;

% Generate channel
[Phi_i,tau_i,v_i] = channel_generation(Fc,vel);
L1 = Q + 1;
L2 = Q + 1 + floor(max(tau_i) / Ts);

% Set up indices for block placement in HDD
u_vec = repmat((0:U-1).',M,1);
m_vec = repelem((0:M-1),1,U).';
indices_base = m_vec*N*U + u_vec;

% Loop through all u,v pairs
HDD = zeros(M*N*U);
for n = 0:N-1
    for k = 0:N-1

        % Generate H matrix block
        HDD_new = gen_Hblock(U,M,Ts,Phi_i,tau_i,v_i,L1,L2,Ambig_Table.vals{n+1,k+1},Ambig_Table.t_range,Ambig_Table.f_range,t_offset);

        % M -> N -> U index assignment
        indices1 = indices_base + k*U + 1;
        indices2 = indices_base + n*U + 1;

        % Add block to stack
        HDD(indices1,indices2) = HDD_new;

    end
end

end

function H_block = gen_Hblock(U,M,Ts,Phi_i,tau_i,v_i,L1,L2,ambig_vals,ambig_t_range,ambig_f_range,t_offset)
% This generates the delay-frequency channel matrix for an ODDM system
% INPUTS:
%   U:              number of frequency taps
%   M:              number of subcarriers
%   Ts:             Symbol period
%   Phi_i:          Multipath coefficients
%   tau_i:          Delay taps
%   v_i:            Doppler taps
%   L1:             Lower M-block range
%   L2:             Upper M-block range
%   ambig_vals:     prerendered cross-ambiguity table values
%   ambig_t_range:  prerendered cross-ambiguity table delay range
%   ambig_f_range:  prerendered cross-ambiguity table frequency range
%
% Coded by Jeremiah Rhys Wimer, 5/19/2025

% Parameters
Fs = 1 / Ts;

% Set up L tap range
range1 = -M + (1:L2);
range2 = -L1:L2;
range3 = M - (L1:-1:1);
L_range = [range1, range2, range3].';

% Find all combinations of m,n,l,k - then only select valid ones
m_range = 0:(M-1);
u_range = 0:(U-1);
l_range = 0:(M-1);
v_range = 0:(U-1);
tap_combos = combvec(m_range,u_range,l_range,v_range).';
l_minus_m = tap_combos(:,3) - tap_combos(:,1);
needed_indices = ismember(l_minus_m,L_range);
needed_combos = tap_combos(needed_indices,:);
m = needed_combos(:,1);
u = needed_combos(:,2);
l = needed_combos(:,3);
v = needed_combos(:,4);

% Find indices in ambiguity values corresponding with the coefficient eq.
ambig_t_locs = (l-m).*Ts - tau_i + t_offset;
ambig_f_locs = (v-u).*Fs - v_i;

tIdxInterp = griddedInterpolant(ambig_t_range   , 1:numel(ambig_t_range), ...
                                'nearest','nearest');
fIdxInterp = griddedInterpolant(ambig_f_range   , 1:numel(ambig_f_range), ...
                                'nearest','nearest');
ambig_t_indices = reshape(tIdxInterp(ambig_t_locs), size(ambig_t_locs));
ambig_f_indices = reshape(fIdxInterp(ambig_f_locs), size(ambig_f_locs));

% Clamp values to stay within bounds
ambig_t_indices = min(max(ambig_t_indices, 1), length(ambig_t_range));
ambig_f_indices = min(max(ambig_f_indices, 1), length(ambig_f_range));

% Convert subscripts to linear indices and select ambiguity values
linear_indices = sub2ind(size(ambig_vals), ambig_t_indices, ambig_f_indices);
ambig_inst = ambig_vals(linear_indices);

% Define all elements of h summation, then sum across 2nd dimension
h_sum = Phi_i .* exp(1j.*2.*pi.*(v_i + u.*Fs).*(l.*Ts-tau_i+t_offset)) .* ambig_inst;
h = sum(h_sum,2);

% Find all linear indices for H matrix
H_block = zeros(M*U);
H_index1 = l*U + v+1;
H_index2 = m*U + u+1;
linear_indices = sub2ind(size(H_block), H_index1, H_index2);

% Assign values from h to the corresponding locations in H
H_block(linear_indices) = h;

end