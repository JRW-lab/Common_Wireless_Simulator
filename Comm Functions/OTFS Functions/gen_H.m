function H = gen_H(T,N,M,Lp,Ln,chn_g,chn_tau,chn_v,shape,alpha,t_offset)
% This generates the time-domain channel matrix for an OTFS system
% INPUTS:
%   Fc: carrier frequency
%   v:  maximum vehicle velocity
%   df: subcarrier spacing
%   N:  number of time symbols
%   M:  number of subcarriers
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024

% rng('default');

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
% q = 4;

% Internal settings
Ts = T / M;

% Check to make sure channel length is still valid
if M <= min(Lp - Ln + 1)
    error("M MUST BE GREATER THAN TOTAL CHANNEL LENGTH")
end

% % METHOD 1 - USES FOR LOOPS
% % Create channel coefficients
% h = zeros(M*N,Lp-Ln+1);
% for k = 0:(N*M)-1
%     count = 0;
%     for l = Ln:Lp
%         count = count + 1;
%         % Render a vector of exponential values
%         exp_vals = exp(1j*2*pi*chn_v*(k-l)*Ts);
% 
%         % Render a vector of cross-ambiguity function values
%         ambig_t_range = l*Ts - chn_tau;
%         ambig_f_range = chn_v;
%         ambig_vals = ambig(ambig_t_range, ambig_f_range, Ts, shape, alpha);
% 
%         % Create sum vector
%         sum_vals = chn_g .* exp_vals .* ambig_vals;
% 
%         % Sum all elements to make value of h
%         h(k+1,count) = sum(sum_vals,"all");
%     end
% end

% % METHOD 2 - USES 2D SPACE
% % Define vector of l given L- and L+
% l = (Ln:Lp).';
% 
% % Premake all possible t and f ranges for h
% ambig_t_range = l*Ts - chn_tau;
% ambig_f_range = ones(Lp-Ln+1,1) .* chn_v;
% 
% % Create channel coefficients
% h = zeros(M*N,Lp-Ln+1);
% for k = 0:(N*M)-1
%     % Render a vector of exponential values
%     exp_vals = exp(1j.*2.*pi.*chn_v.*(k-l).*Ts);
% 
%     % Render a vector of cross-ambiguity function values
%     ambig_vals = ambig(ambig_t_range, ambig_f_range, Ts, shape, alpha);
% 
%     % Create sum vector
%     sum_vals = chn_g .* exp_vals .* ambig_vals;
% 
%     % Sum all elements to make value of h
%     h(k+1,:) = sum(sum_vals,2).';
% end

% METHOD 3 - USES 3D SPACE
% Define vector of l given L- and L+
% Ln_h_gen = min(Ln);
% Lp_h_gen = max(Lp);
l = (Ln:Lp).';

% Make all possible t and f ranges for h
ambig_t_range = (l*Ts - chn_tau + t_offset) .* ones(Lp-Ln+1,length(chn_g),N*M);
ambig_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),N*M);

% Make all possible exponential values for h
k = reshape(0:(N*M-1), [1, 1, N*M]);
exp_vals = exp(1j.*2.*pi.*chn_v.*((k-l).*Ts + t_offset));

% Make ambiguity values and sum to make h
ambig_vals = ambig(ambig_t_range, ambig_f_range, Ts, shape, alpha);
sum_vals = chn_g .* exp_vals .* ambig_vals;
h = squeeze(sum(sum_vals,2)).';

% Create channel matrix from coefficients
H = zeros(M*N);
H(:,1:Lp-Ln+1) = fliplr(h(:,1:Lp-Ln+1));
for k = 1:M*N
    H(k,:) = circshift(H(k,:),k-Lp-1);
end
