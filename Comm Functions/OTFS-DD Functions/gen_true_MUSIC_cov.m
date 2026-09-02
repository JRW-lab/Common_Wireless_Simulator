function R_x = gen_true_MUSIC_cov(num_pilots,pilot_spacing,Ts,l,chn_tau,chn_v,N0,ambig_vals,ambig_t_range,ambig_f_range)

% Parameters
pdp = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];        % in dB
pow_prof = 10.^(pdp/10);                                    % dB to watt
sigma2_p = pow_prof/sum(pow_prof);                          % normalization

% Set up range of (k1-k2) and exponential part
k_diff_range = (((num_pilots-1)*pilot_spacing):-pilot_spacing:-((num_pilots-1)*pilot_spacing)).';
omega = 2*pi*chn_v*Ts;
exp_mat = exp(1j.*omega.*k_diff_range);

% Find absolute ambiguity values
t_vals = l*Ts - chn_tau;
f_vals = chn_v;
t_val_diff = abs(ambig_t_range.' - t_vals).^2;
[~,t_idx] = min(t_val_diff,[],1);
f_val_diff = abs(ambig_f_range.' - f_vals).^2;
[~,f_idx] = min(f_val_diff,[],1);
ambig_vals_sel = abs(diag(ambig_vals(t_idx,f_idx)).').^2;

% Create sigma2_p_tilde
sigma2_p_tilde = ambig_vals_sel .* sigma2_p;

% Create the extended row of R_x
cov_whole_row = sum(exp_mat .* sigma2_p_tilde,2).';

% Build R_x
R_x = zeros(num_pilots);
for tap_idx = 0:num_pilots-1
    shifted_row_ext = circshift(cov_whole_row,num_pilots+tap_idx);
    R_x(tap_idx+1,:) = shifted_row_ext(1:num_pilots);
end

% Add noise component along main diagonal
R_x = R_x + N0 * eye(num_pilots);
