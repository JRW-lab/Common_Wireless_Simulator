function metrics_aux = merge_metrics_aux(old_aux, new_aux)
% MERGE_METRICS_AUX Combine running aggregates for the stability check.
%
%   metrics_aux = merge_metrics_aux(old_aux, new_aux)
%
%   old_aux and new_aux are structs with running aggregate fields.
%   Binomial-style counters (total_bit_errors, total_bits, ...) are added
%   directly. Continuous statistics tracked via Welford's algorithm
%   (mean/M2/count) are merged using the parallel-batch update formula.
%
%   If old_aux is empty (first run), new_aux is returned unchanged.

% Handle first run
if isempty(old_aux)
    metrics_aux = new_aux;
    return;
end

metrics_aux = struct();

% Binomial-style additive counters
add_fields = {'total_bit_errors', 'total_bits', ...
              'total_sym_errors', 'total_symbols', ...
              'total_frm_errors', 'total_frames'};
for i = 1:numel(add_fields)
    f = add_fields{i};
    if isfield(new_aux, f)
        old_val = 0;
        if isfield(old_aux, f)
            old_val = old_aux.(f);
        end
        metrics_aux.(f) = old_val + new_aux.(f);
    elseif isfield(old_aux, f)
        metrics_aux.(f) = old_aux.(f);
    end
end

% Continuous statistics via Welford's parallel batch merge.
% Merge every continuous field block named <name>_count / <name>_mean / <name>_M2,
% detected automatically from new_aux (generic; covers t_RXfull, recon_mse, ...).
add_set = { 'total_bit_errors', 'total_bits', ...
            'total_sym_errors', 'total_symbols', ...
            'total_frm_errors', 'total_frames' };
cont_names = {};
new_fields = fieldnames(new_aux);
for i = 1:numel(new_fields)
    f = new_fields{i};
    if endsWith(f, '_count') && ~ismember(f, add_set)
        name = extractBefore(f, '_count');
        if isfield(new_aux, [name '_mean']) && isfield(new_aux, [name '_M2'])
            cont_names{end+1} = name; %#ok<AGROW>
        end
    end
end

for c = 1:numel(cont_names)
    name = cont_names{c};
    % New values
    n_new = new_aux.([name '_count']);
    mean_new = new_aux.([name '_mean']);
    M2_new = new_aux.([name '_M2']);
    % Old values (default to zero state if absent)
    if isfield(old_aux, [name '_count'])
        n_old = old_aux.([name '_count']);
        mean_old = old_aux.([name '_mean']);
        M2_old = old_aux.([name '_M2']);
    else
        n_old = 0;
        mean_old = 0;
        M2_old = 0;
    end
    n_total = n_old + n_new;
    if n_total > 0
        delta = mean_new - mean_old;
        metrics_aux.([name '_count']) = n_total;
        metrics_aux.([name '_mean']) = (n_old * mean_old + n_new * mean_new) / n_total;
        metrics_aux.([name '_M2']) = M2_old + M2_new + delta^2 * n_old * n_new / n_total;
    else
        metrics_aux.([name '_count']) = 0;
        metrics_aux.([name '_mean']) = 0;
        metrics_aux.([name '_M2']) = 0;
    end
end

end
