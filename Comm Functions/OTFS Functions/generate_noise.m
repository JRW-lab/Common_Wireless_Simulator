function result = generate_noise(noise_covariance_matrix)

[U_z, D_z, ~] = svd(noise_covariance_matrix);
R_z_half = (U_z * sqrt(D_z));

syms_per_f = size(R_z_half,1);

n = sqrt(1 / 2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
result = R_z_half * n;

end