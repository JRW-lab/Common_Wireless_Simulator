function [TX_bit,TX_sym,s] = gen_data(bin_order,S,num_data_syms,rs)
% Generate random bits, symbol equivalents, and modulated symbols.
%     bin_order = binary order in the form of [0,0]
%                                             [0,1]
%                                             [1,1]
%                                             [1,0]
%             S = alphabet for symbols in the form a vector, where row n of bin_order
%                 corresponds with a symbol of index n
% num_data_syms = number of desired symbols
%            rs = optional RandStream for the symbol draw (Phase 3); global stream when absent
%
% The outputs are in the form of bits, symbol indices, and complex-value
% symbols.


% Create TX symbols
if nargin < 4 || isempty(rs)
    TX_sym = randi(length(S),num_data_syms,1)-1;
else
    TX_sym = randi(rs,length(S),num_data_syms,1)-1;
end

% Create TX bits
TX_bit = bin_order(TX_sym+1,:);

% Create modulated symbols s
s = S(TX_sym+1);