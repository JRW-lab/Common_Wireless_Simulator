function [parameters, paramHash, is_real] = build_point_params(prof, primary_val, config_idx)
%BUILD_POINT_PARAMS  The ONE definition of "what parameter struct does this
%test point have, and therefore what is its param_hash".
%
%   [parameters, paramHash, is_real] = build_point_params(prof, primary_val, config_idx)
%
%   prof         a profile struct from saved_profiles()
%   primary_val  a value of prof.primary_var
%   config_idx   index into prof.configs
%   is_real      false when this (primary_val, config) pair is NOT a test
%                point -- i.e. a flat config away from its anchor. Callers
%                MUST skip those rather than collect or query them.
%
% WHY THIS FUNCTION EXISTS. Reconstructing this merge by hand is the single
% most repeated bug in this project: Bug Log #6a and #12 are both "a
% hand-built param struct hashed differently from what sim_head wrote, so
% the lookup silently returned no data and collection forked into an
% orphaned row". Nothing errors when that happens. Before this function the
% same ~15 lines were duplicated in sim_head, collect_profile, wait_profile,
% render_profile, analyze_ngs, render_figure5 and the GUI's Frames button --
% every one of them a place for the merge order, the ODDM `U` scrub or the
% flat-config rules to drift out of step.
%
% NOT a caller of this function, deliberately: collectors/hash_regression.m
% keeps its own independent copy of the enumeration, so that a bug
% introduced HERE is caught by its before/after diff instead of being
% mirrored by it. Do not "tidy" that duplication away.
%
% FLAT CONFIGS. A config listed in prof.flat_configs holds one measurement
% that does not depend on the primary variable (see sim_head.m's FLAT
% CONFIGS note). Two modes:
%   - without prof.flat_defaults: built like any other point but only at
%     prof.flat_anchor, which therefore enters the hash and is required.
%   - with prof.flat_defaults: built from that struct instead, and the
%     primary variable is NOT injected, so the swept parameter does not
%     appear in the hash at all. This is what lets a flat config match a
%     measurement collected by a separate standalone profile.

flat_cols = [];
if isfield(prof,'flat_configs') && ~isempty(prof.flat_configs)
    flat_cols = unique(prof.flat_configs(:)).';
end
is_flat = ismember(config_idx, flat_cols);
has_flat_defaults = isfield(prof,'flat_defaults') && ~isempty(prof.flat_defaults);

% Is this pair a real point?
is_real = true;
if is_flat
    if isfield(prof,'flat_anchor') && ~isempty(prof.flat_anchor)
        is_real = (primary_val == prof.flat_anchor);
    elseif has_flat_defaults
        is_real = (primary_val == prof.primary_vals(1));
    else
        error("build_point_params:missingFlatAnchor", ...
            "profile declares flat_configs but neither flat_anchor nor flat_defaults.");
    end
end
if ~is_real
    parameters = [];
    paramHash = '';
    return;
end

% ---- sim_head.m's merge order, and nothing else ----
if is_flat && has_flat_defaults
    parameters = prof.flat_defaults;
    inject_primary = false;
else
    parameters = prof.default_parameters;
    inject_primary = true;
end
if isfield(prof,'MUSIC_settings')
    parameters = mergestructs(parameters, prof.MUSIC_settings);
end
if isfield(prof,'vehicle_motion_settings')
    parameters = mergestructs(parameters, prof.vehicle_motion_settings);
end
if inject_primary
    parameters.(prof.primary_var) = primary_val;
end
cfg = prof.configs{config_idx};
fn = fieldnames(cfg);
for k = 1:numel(fn)
    parameters.(fn{k}) = cfg.(fn{k});
end

% ---- per-system field scrubbing, identical to sim_head.m ----
switch parameters.system_name
    case {"ODDM","OTFS"}
        parameters = rmfield(parameters,'U');
    case "OFDM"
        parameters = rmfield(parameters,'N');
        parameters = rmfield(parameters,'U');
        parameters = rmfield(parameters,'shape');
        parameters = rmfield(parameters,'alpha');
        parameters = rmfield(parameters,'Q');
end
% sim_head.m has a further block here guarded by
%     if exist("parameters.shape",'var')
% which is ALWAYS FALSE -- exist() takes a variable name, not a field
% expression, so "parameters.shape" never resolves and the body never runs.
% It is therefore dead code, and reproducing it faithfully means omitting
% it. Do NOT "fix" that guard in sim_head: the body removes `alpha` for
% non-rrc pulses and forces Q=1 for rect/ideal, so making it live would
% change the param_hash of every affected point and orphan its data. Logged
% rather than silently inherited.

% MUSIC profiles carry a large additional clean-up block in sim_head that is
% deliberately NOT duplicated here. Rather than guess and produce a hash
% that merely looks plausible -- the exact failure this function exists to
% prevent -- refuse the case outright.
if isfield(parameters,'channel_estimation_method')
    error("build_point_params:unsupportedMusicProfile", ...
        "This profile uses channel_estimation_method (a MUSIC-style profile). " + ...
        "sim_head applies further field clean-up for those that is not " + ...
        "replicated here, so the hash would be wrong. Extend this function " + ...
        "to match sim_head before using it on such a profile.");
end

[~, paramHash] = jsonencode_sorted(parameters);
end
