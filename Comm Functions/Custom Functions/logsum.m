function y = logsum(x)
%LOGSUM Stable computation of log(sum(exp(x))) returning scalar.
%
% x is a vector in the log-domain
% y = log(exp(x1)+exp(x2)+...+exp(xN))

    if isempty(x)
        error('Input cannot be empty.');
    end

    % Use the max-trick for numerical stability
    m = max(x);
    y = m + log(sum(exp(x - m)));
end
