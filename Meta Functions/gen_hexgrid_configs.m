function configs = gen_hexgrid_configs(N_sym, M_min, N_min, U_min)

configs = {};  % Initialize cell array to store structs
idx = 1;       % Index for configs

% Loop over all divisors of N_sym for M
for M = M_min:N_sym
    if mod(N_sym, M) ~= 0
        continue;
    end
    rem1 = N_sym / M;

    % Loop over all divisors of the remainder for N
    for N = N_min:rem1
        if mod(rem1, N) ~= 0
            continue;
        end
        U = rem1 / N;

        % Check U constraint
        if U >= U_min && mod(U,1) == 0
            configs{idx} = struct('M', M, 'N', N, 'U', U);
            idx = idx + 1;
        end
    end
end
end