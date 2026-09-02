function [H,h] = gen_H(Ts,syms_per_f,Lp,Ln,chn_g,chn_tau,chn_v,shape,alpha,t_offset)
% This generates the time-domain channel matrix for an OTFS system
% INPUTS:
%   Fc: carrier frequency
%   v:  maximum vehicle velocity
%   df: subcarrier spacing
%   N:  number of time symbols
%   M:  number of subcarriers
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024 - modified 1/30/2026

% METHOD 3 - USES 3D SPACE
l = (Ln:Lp).';

% Make all possible t and f ranges for h
ambig_t_range = (l*Ts - chn_tau + t_offset) .* ones(Lp-Ln+1,length(chn_g),syms_per_f);
ambig_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),syms_per_f);

% Make all possible exponential values for h
k = reshape(0:(syms_per_f-1), [1, 1, syms_per_f]);
exp_vals = exp(1j.*2.*pi.*chn_v.*((k-l).*Ts + t_offset));

% Make ambiguity values and sum to make h
ambig_vals = ambig(ambig_t_range, ambig_f_range, Ts, shape, alpha);
sum_vals = chn_g .* exp_vals .* ambig_vals;
h = squeeze(sum(sum_vals,2)).';

% Create channel matrix from coefficients
H = zeros(syms_per_f);
H(:,1:Lp-Ln+1) = fliplr(h(:,1:Lp-Ln+1));
for k = 1:syms_per_f
    H(k,:) = circshift(H(k,:),k-Lp-1);
end
