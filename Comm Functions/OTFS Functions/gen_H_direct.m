function H = gen_H_direct(T,N,M,Lp,Ln,chn_g,chn_v,ambig_table,tap_t_range,tap_f_range,t_offset)
% This generates the time-domain channel matrix for an OTFS system
% INPUTS:
%   Fc: carrier frequency
%   v:  maximum vehicle velocity
%   df: subcarrier spacing
%   N:  number of time symbols
%   M:  number of subcarriers
%
% Coded by Jeremiah Rhys Wimer, 3/24/2024

% Internal settings
Ts = T / M;

% Check to make sure channel length is still valid
if M <= min(Lp - Ln + 1)
    error("M MUST BE GREATER THAN TOTAL CHANNEL LENGTH")
end

% METHOD 3 - USES 3D SPACE
% Define vector of l given L- and L+
l = (Ln:Lp).';

% Make all possible exponential values for h
k = reshape(0:(N*M-1), [1, 1, N*M]);
exp_vals = exp(1j.*2.*pi.*chn_v.*((k-l).*Ts + t_offset));

% Make ambiguity values and sum to make h
ambig_vals = 0 * tap_t_range;
for k = 1:size(tap_t_range,1)
    for l = 1:size(tap_t_range,2)
        for m = 1:size(tap_t_range,3)
            try
                ambig_vals(k,l,m) = ambig_table(tap_t_range(k,l,m),tap_f_range(k,l,m));
            catch
                error("AMBIG TABLE ERROR");
            end
        end
    end
end
sum_vals = chn_g .* exp_vals .* ambig_vals;
h = squeeze(sum(sum_vals,2)).';

% Create channel matrix from coefficients
H = zeros(M*N);
H(:,1:Lp-Ln+1) = fliplr(h(:,1:Lp-Ln+1));
for k = 1:M*N
    H(k,:) = circshift(H(k,:),k-Lp-1);
end
