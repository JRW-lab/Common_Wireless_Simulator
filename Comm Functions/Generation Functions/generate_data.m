function [TXdata,s] = generate_data(S,num_data_bits)

M = length(S);

TXdata = randi(M, 1, num_data_bits)';
s = S(TXdata);

end