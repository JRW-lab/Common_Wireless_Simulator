function bin_vec = data2bin(data,M)

% % Make only one bit difference between numbers
% if M == 4
%     data_copy = data;
%     for k = 1:length(data)
%         if data(k) == 2
%             data_copy(k) = 3;
%         elseif data(k) == 3
%             data_copy(k) = 2;
%         end
%     end
%     data = data_copy;
% end

bin_mat = dec2bin(data);
[N,L] = size(bin_mat);

bin_vec = repmat(' ', 1, N*L);
count = 0;
for n = 1:N
    for l = 1:L
        count = count + 1;
        bin_vec(count) = bin_mat(n,l);
    end
end

bin_vec = bin_vec.';

end