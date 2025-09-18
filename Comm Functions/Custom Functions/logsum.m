function y = logsum(x)
% function y = logsum(x)
% 
% x is a vector in the log-domain
% y = log(exp(x1)+exp(x2)+...+exp(xN))
%
% algorithm:
% log(exp(x1) + exp(x2)) = log(exp(x1))+log(1+exp(x2-x1)) = x1 + log[1+exp(x2-x1)]
% the results are calculated recursively
%


if length(x) == 1
	y = x;
elseif length(x) > 2
	temp_y = logsum(x(1:2));
	y = logsum([temp_y, x(3:end)]);
elseif length(x) == 2
	y = max(x) + log(1+exp(min(x)-max(x)));
else
	error('something is wrong!');
end
 


