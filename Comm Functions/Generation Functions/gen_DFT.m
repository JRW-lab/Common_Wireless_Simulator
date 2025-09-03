function F_N = gen_DFT(size)
% This is a function for generating the normalized N-point Discrete Fourier
% Transform matrix
omega = exp(-1j * 2*pi / (size));
F_N = zeros(size);
for m = 1:1:size
    for n = 1:1:size
        F_N(m,n) = omega^((m-1) * (n-1));
    end
end

F_N = F_N / sqrt(size);

end