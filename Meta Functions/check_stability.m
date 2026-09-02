function [is_stable, current_value, ci_width] = check_stability(metrics_aux, data_type, relative_tolerance, min_frames, confidence)
% CHECK_STABILITY Determine if a data point has converged.
%
%   [is_stable, current_value, ci_width] = check_stability(metrics_aux, ...
%       data_type, relative_tolerance, min_frames, confidence)
%
%   metrics_aux is a struct of running aggregates (see merge_metrics_aux).
%   For binomial metrics (BER/SER/FER) a Wilson score interval is used; for
%   continuous metrics a t-distribution interval is built from Welford's
%   running mean/M2/count.
%
%   A point is declared stable when its confidence-interval half-width is
%   smaller than relative_tolerance * |current_value| (dynamic threshold
%   that shrinks as the metric's decade shrinks).

is_stable = false;
current_value = NaN;
ci_width = Inf;

if isempty(metrics_aux) || ~isfield(metrics_aux, 'total_frames')
    return;
end

N_frames = metrics_aux.total_frames;
if N_frames < min_frames
    return;
end

z = norminv(1 - (1 - confidence) / 2);

switch data_type
    case "BER"
        if ~isfield(metrics_aux, 'total_bits') || ~isfield(metrics_aux, 'total_bit_errors')
            return;
        end
        total_trials = metrics_aux.total_bits;
        total_errors = metrics_aux.total_bit_errors;
        if total_trials == 0
            return;
        end
        p_hat = total_errors / total_trials;
        ci_width = z * sqrt(p_hat * (1 - p_hat) / total_trials) ...
            / (1 + z^2 / total_trials);
        current_value = p_hat;

    case "SER"
        if ~isfield(metrics_aux, 'total_symbols') || ~isfield(metrics_aux, 'total_sym_errors')
            return;
        end
        total_trials = metrics_aux.total_symbols;
        total_errors = metrics_aux.total_sym_errors;
        if total_trials == 0
            return;
        end
        p_hat = total_errors / total_trials;
        ci_width = z * sqrt(p_hat * (1 - p_hat) / total_trials) ...
            / (1 + z^2 / total_trials);
        current_value = p_hat;

    case "FER"
        if ~isfield(metrics_aux, 'total_frames') || ~isfield(metrics_aux, 'total_frm_errors')
            return;
        end
        total_trials = metrics_aux.total_frames;
        total_errors = metrics_aux.total_frm_errors;
        if total_trials == 0
            return;
        end
        p_hat = total_errors / total_trials;
        ci_width = z * sqrt(p_hat * (1 - p_hat) / total_trials) ...
            / (1 + z^2 / total_trials);
        current_value = p_hat;

    otherwise
        % Continuous metric via Welford's running statistics
        name = char(data_type); % e.g. "t_RXfull"
        if ~isfield(metrics_aux, [name '_count']) || ~isfield(metrics_aux, [name '_mean']) ...
                || ~isfield(metrics_aux, [name '_M2'])
            return;
        end
        n = metrics_aux.([name '_count']);
        if n < 2
            return;
        end
        mu = metrics_aux.([name '_mean']);
        M2 = metrics_aux.([name '_M2']);
        current_value = mu;
        std_val = sqrt(M2 / (n - 1));
        se = std_val / sqrt(n);
        t_crit = tinv(1 - (1 - confidence) / 2, n - 1);
        ci_width = t_crit * se;
end

% Dynamic threshold shrinks with the metric's decade
if current_value > 0
    threshold = current_value * relative_tolerance;
else
    threshold = eps;
end
is_stable = (ci_width < threshold);
end
