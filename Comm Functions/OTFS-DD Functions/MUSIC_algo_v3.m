function [R_x,s_est,v_est,recon_norm,P_hat] = MUSIC_algo_v3(sys_data,X,MUSIC_settings,debug_settings)

% Import debug settings
use_true_x_cov = debug_settings.use_true_x_cov;
tap_idx = debug_settings.tap_idx;
chn_tau = debug_settings.chn_tau;
chn_v = debug_settings.chn_v;
Ln = debug_settings.Ln;
N0 = debug_settings.N0;
ambig_vals = debug_settings.ambig_vals;
ambig_t_range = debug_settings.ambig_t_range;
ambig_f_range = debug_settings.ambig_f_range;

% Import MUSIC settings
enforce_toeplitz = MUSIC_settings.enforce_toeplitz;
shrinkage_alpha = MUSIC_settings.shrinkage_alpha;
average_cov = MUSIC_settings.average_cov;
cov_epsilon = MUSIC_settings.cov_epsilon;
R_x_last = MUSIC_settings.R_x_last;

% Import system data
Ts = sys_data.Ts;
num_pilots = sys_data.num_pilots;
pilot_idx = sys_data.pilot_idx;
frame_split_divisor = sys_data.frame_split_divisor;
vel = sys_data.vel;
Fc = sys_data.Fc;
v_res = sys_data.v_res;

% Set variables for MUSIC
pilot_spacing = unique(diff(pilot_idx));
split_size = num_pilots/frame_split_divisor;
v_max = (vel * (1000/3600)*Fc) / physconst('LightSpeed');
v_max_rounded = 100 * ceil(v_max / 100);
v_range = -v_max_rounded:v_res:v_max_rounded;
power_vec = pilot_idx-1;
e = exp(1j*2*pi*Ts.*power_vec(1:frame_split_divisor:end).*v_range);

% Generate covariance and select noise subspace
if ~use_true_x_cov

    % Perform frame splitting
    R_x = zeros(split_size);
    for i = 1:frame_split_divisor
        X_sec = X(i:frame_split_divisor:end,:);
        R_x_new = X_sec * X_sec' / size(X_sec,2);
        R_x = R_x + R_x_new;
    end
    R_x = R_x / frame_split_divisor;

    % Perform shrinkage
    R_x = (1 - shrinkage_alpha) * R_x + shrinkage_alpha * (trace(R_x) / split_size) * eye(split_size);

    % Enforce toeplitz structure
    if enforce_toeplitz
        R_toeplitz = zeros(split_size);
        for d = -(split_size-1):(split_size-1)
            vals = diag(R_x,d);
            avg = mean(vals);
            R_toeplitz = R_toeplitz + diag(avg*ones(length(vals),1), d);
        end
        R_x = R_toeplitz;
    end

else
    R_x = gen_true_MUSIC_cov(split_size,pilot_spacing*frame_split_divisor,Ts,tap_idx-1+Ln,chn_tau,chn_v,N0,ambig_vals,ambig_t_range,ambig_f_range);
end

% Average with last covariance matrix
if average_cov
    R_x = (1-cov_epsilon) * R_x + cov_epsilon * R_x_last;
end

% Perform eigenvalue decomposition on R_x
[U,D] = eig(R_x);
[~, ind] = sort(diag(D), 'descend');
Us = U(:,ind);

% Test reconstruction quality of all possible number of paths
num_paths = 10; % More paths than this seems unlikely
recon_norm_mat = zeros(num_paths,1);
v_est_cell = cell(num_paths,1);
P_hat_cell = cell(num_paths,1);
A_est_cell = cell(num_paths,1);
s_est_cell = cell(num_paths,1);
for p_est = 1:num_paths

    % Select noise subspace
    Un = Us(:,(p_est+1):end);

    % Sweep through all possible Doppler values and generate costs
    norm_mat = Un' * e;
    d_sqrd = sum(abs(norm_mat).^2,1);
    P_hat_cell{p_est} = 1 ./ d_sqrd;

    % Estimate doppler shifts
    [~,v_est_idx] = findpeaks(real(P_hat_cell{p_est}),'NPeaks',p_est,'SortStr','descend');
    v_est_cell{p_est} = v_range(v_est_idx).';

    % Find reconstruction difference norm
    [recon_norm_mat(p_est),A_est_cell{p_est},s_est_cell{p_est}] = calc_recon_norm(X,Ts,power_vec(1:num_pilots),v_est_cell{p_est});

end

% Select best fitting Doppler set
[~,idx_sel] = min(recon_norm_mat);
v_est = v_est_cell{idx_sel};
recon_norm = recon_norm_mat(idx_sel);
P_hat = P_hat_cell{idx_sel};
s_est = s_est_cell{idx_sel};
