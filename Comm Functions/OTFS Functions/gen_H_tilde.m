function H_tilde = gen_H_tilde(H,M,N,Lp,Ln,F_N)
% This takes a time domain channel matrix for OTFS and permutes it for
% easier equalization. This can be achieved through block transformations
% that are more efficient for MATLAB.
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024

% % Debugging inputs
% clear; clc;
% obj = comms_obj_OTFS;
% Fc = obj.Fc;
% v = obj.v_vel;
% df = obj.sbcar_spacing;
% N = obj.N_tsyms;
% M = obj.M_sbcars;
% shape = obj.filter;
% alpha = obj.rolloff;
% q = obj.q;
% [H,Ln,Lp] = gen_H(Fc,v,df,N,M,shape,alpha,q);
% F_N = gen_DFT(N);

% rng('default');

% % METHOD 1
% tic;
% % Define range that l can take (0:L)
% l_vec = Ln:Lp;
% 
% % Render H_tilde each block at a time
% H_tilde = zeros(N*M);
% for u = 0:M-1
%     for v = 0:M-1
%         % Check logic statement from paper
%         if ismember(v,mod(u-l_vec,M))
%             % Create block from given H matrix
%             H_block = zeros(N);
%             for m = 0:N-1
%                 for n = 0:N-1
%                     H_block(m+1,n+1) = H((m*M)+u+1,(n*M)+v+1);
%                 end
%             end
% 
%             % Take DFT of newly created block
%             H_tilde_block = F_N * H_block * F_N';
% 
%             % Put block in corresponding slot of H_tilde
%             H_tilde((u*N)+1:(u+1)*N,(v*N)+1:(v+1)*N) = H_tilde_block;
%         end
%     end
% end
% toc
% % end

% METHOD 2
% Define range that l can take (Ln:Lp)
l_vec = reshape(Ln:Lp,[1,1,Lp-Ln+1]);

% Loop over u and v together using meshgrid
[v_grid, u_grid] = meshgrid(0:M-1);

% Compute mod(u - l_vec, M) once outside the loop
mod_values = mod(u_grid - l_vec, M);

% Find the indices where v is a member of the mod values
logic_map = not(sum(v_grid == mod_values,3) > 0.5);
u_grid(logic_map) = NaN;
v_grid(logic_map) = NaN;
u_grid = rmmissing(u_grid(:));
v_grid = rmmissing(v_grid(:));

% Render H_tilde each block at a time
H_tilde = zeros(N*M);
for k = 1:numel(u_grid)
    % Set current indices
    u = u_grid(k);
    v = v_grid(k);

    % Create block from given H matrix using indexing
    H_block = H(u+1:M:N*M, v+1:M:N*M);

    % Take DFT of newly created block
    H_tilde_block = F_N * H_block * F_N';

    % Put block in corresponding slot of H_tilde
    H_tilde((u*N)+1:(u+1)*N, (v*N)+1:(v+1)*N) = H_tilde_block;
end
