function [TXdata,s] = generate_data(S,num_data_syms)

M = length(S);

TXdata = randi(M, 1, num_data_syms)';
s = S(TXdata);

end