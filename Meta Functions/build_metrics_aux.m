function new_aux = build_metrics_aux(new_frames, frame_data)
% BUILD_METRICS_AUX Build running aggregate metrics_aux from frame_data.
%
%   new_aux = build_metrics_aux(new_frames, frame_data)
%
%   frame_data is the per-frame struct returned by a sim_fun_* (bit_errors,
%   sym_errors, frm_errors, t_RXfull, bits_per_frame, syms_per_frame, and any
%   optional per-frame continuous fields such as recon_mse).
%
%   The returned new_aux contains binomial additive counters plus Welford
%   running statistics (_count/_mean/_M2) for every continuous per-frame
%   field present, so check_stability can key off any of them.
%
%   Continuous-field validity follows the established convention: a per-frame
%   value of 0 means "not computed on this frame" (e.g. t_RXfull on a skipped
%   frame, or recon_mse on an all-pilot frame) and is excluded from the
%   running mean/M2 but still counted toward the frame total via the additive
%   counters.

% Binomial additive counters
new_aux = struct();
% Use the true number of data frames if the sim_fun reports it (some systems
% interleave pilot-only frames with data frames, e.g. MUSIC OTFS). Defaults
% to new_frames when frame_data.n_data_frames is not supplied.
if isfield(frame_data, 'n_data_frames') && ~isempty(frame_data.n_data_frames)
    n_frames_eff = frame_data.n_data_frames;
else
    n_frames_eff = new_frames;
end
new_aux.total_bit_errors = sum(frame_data.bit_errors);
new_aux.total_bits       = n_frames_eff * frame_data.bits_per_frame;
new_aux.total_sym_errors = sum(frame_data.sym_errors);
new_aux.total_symbols    = n_frames_eff * frame_data.syms_per_frame;
new_aux.total_frm_errors = sum(frame_data.frm_errors);
new_aux.total_frames     = n_frames_eff;

% Continuous Welford stats for every continuous per-frame field present.
cont_fields = {'t_RXfull', 'recon_mse'};
for i = 1:numel(cont_fields)
    name = cont_fields{i};
    if isfield(frame_data, name)
        v = double(frame_data.(name));
        if size(v,1) == 1 && size(v,2) > 1
            v = v(:);
        end
        valid = (v > 0) & ~isnan(v);
        n = sum(valid);
        if n > 0
            vv = v(valid);
            new_aux.([name '_count']) = n;
            new_aux.([name '_mean'])  = mean(vv);
            new_aux.([name '_M2'])    = sum((vv - mean(vv)).^2);
        else
            new_aux.([name '_count']) = 0;
            new_aux.([name '_mean'])  = 0;
            new_aux.([name '_M2'])    = 0;
        end
    end
end

end
