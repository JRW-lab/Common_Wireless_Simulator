function [x_hat,iters,t_RXiter,t_RXfull] = equalizer_MMSE(y,H,Es,N0)

% Start runtime
tStartRXiter = tic;

% Create MMSE matrix and the corresponding output
A = H' * H + (N0/Es) * eye(length(y));
L = chol(A,"lower");
b = H' * y;
x_hat = L'\(L\b);

% Stop runtime
t_RXiter = toc(tStartRXiter);
t_RXfull = t_RXiter;
iters = 1;