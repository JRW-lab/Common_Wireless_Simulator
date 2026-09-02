function [recon_norm,A_est,s_est] = calc_recon_norm(X,Ts,power_vec,v_est)

% Recreate A from Doppler estimates to get s vectors
A_est = exp(1j*2*pi*Ts.*power_vec.*(sort(v_est).'));
s_est = lsqminnorm(A_est,X);

% Reconstruct X with this frame's Doppler estimates
X_recon = A_est * s_est;
recon_norm = sum(abs(X_recon - X).^2,"all");
