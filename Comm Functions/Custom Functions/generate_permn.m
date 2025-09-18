function states = generate_permn(M, L)
% FUNCTION NOISE = generate_permn(M, L)
%    M: Modulation Order
%    L: Channel Length
%    Generates all possible combinations of M for L channel taps

% INPUTS - FOR DEBUGGING
% clear; clc;
% M = 2;
% L = 5;

states = ones(M^L,L);
for l = 1:L
    for i = 1:M^L
        states(i,l) = floor(mod((i-1) / (M^(L-l)),M) + 1);
    end
end

end