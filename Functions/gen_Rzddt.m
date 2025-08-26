function [R_half,R_zdd_tilde] = gen_Rzddt(T,N,M,q,res_discretization,shape,alpha)

% T = 1/15000;
% N = 16;
% M = 64;
% q = 2;
% res_discretization = 10;
% shape = "sinc";
if not(exist('alpha','var'))
    alpha = 1;
end

% Create duration of data symbol
Ts = T/M;

% Generate needed matrices
F_N = gen_DFT(N);
F_M = gen_DFT(M);
Gamma_MN = gen_Gamma_MN(M,N);

% Generate all possible noise taps
R_taps = zeros(2*N-1,2*M-1);
N_indices = -N+1:N-1;
M_indices = -M+1:M-1;
for n = 1:length(N_indices)
    for m = 1:length(M_indices)
        R_taps(n,m) = ambig_direct((N_indices(n))*Ts,(M_indices(m))/Ts,Ts,shape,alpha,q,res_discretization);
    end
end

% Assign correct noise taps to covariance matrix
R_ztf = zeros(N*M);
% fprintf("| Progress:")
% for i = 1:M-12
%     fprintf(" ");
% end
% fprintf("|\n");
for m1 = 0:M-1
    for m2 = 0:M-1
        for n1 = 0:N-1
            for n2 = 0:N-1
                n_num = n1-n2;
                m_num = m2-m1;
                n_index = N_indices == n_num;
                m_index = M_indices == m_num;
                R_ztf(m1*N+n1+1,m2*N+n2+1) = R_taps(n_index,m_index);
            end
        end
    end
    % fprintf("x")
end

% Convert covariance matrix into needed domain
R_zdd = kron(F_N',F_M)' * R_ztf * kron(F_N',F_M);
R_zdd_tilde = Gamma_MN' * R_zdd * Gamma_MN;

% Generate "half" matrix of covariance
[U_z, D_z, ~] = svd(R_zdd_tilde);
R_half = U_z * sqrt(D_z);
