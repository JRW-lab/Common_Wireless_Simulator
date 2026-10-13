function flat = flatten_params(parameters)
%FLATTEN_PARAMS  Canonical single-level form of a test point's parameters.
%
%   flat = flatten_params(parameters)
%
%   Returns a struct with NO sub-structures, so that the same physical
%   setting has one name no matter where a profile put it. Sorting and
%   hashing are left to jsonencode_sorted; this function only decides names.
%
%   WHY. param_hash used to be taken over the nested struct, so placement
%   leaked into the key: Profile 9 set csi_settings.Pmax = 9, Profile 10
%   sweeps a top-level Pmax = 9, the sim_fun treats the two identically
%   (sim_fun_ODDM_DDRELAXSI_PTMMSE's top-level hooks), and yet they hashed
%   to different rows. Every such pair is a duplicate datum. Keying on this
%   flat form makes placement irrelevant.
%
%   NAMING RULES -- they mirror the sweep hooks the sim_funs already read,
%   so a swept top-level field and a config sub-field land on one name:
%
%     csi_settings.<f>   ->  <f>         (top-level gamma_p, NGS, PiTau,
%                                         PiNu and Pmax hooks use bare names)
%     ldpc_settings.<f>  ->  ldpc_<f>    (resolve_ldpc_settings' aliases:
%                                         ldpc_rate, ldpc_N, ldpc_bp_iters)
%                            colWeight -> ldpc_colweight, matching that
%                                         function's alias exactly
%
%   ldpc_settings being PRESENT turns LDPC on (resolve_ldpc_settings), so a
%   present struct always yields ldpc_enable, defaulting to true. Without
%   this an empty ldpc_settings would flatten to nothing and collide with
%   an uncoded run.
%
%   ERRORS rather than guesses, because a wrong guess here produces a
%   plausible hash for the wrong point and nothing downstream notices:
%     * a sub-structure with no naming rule (add one below);
%     * two sources for one flat name (e.g. top-level Pmax AND
%       csi_settings.Pmax), which the sim_funs reject as ambiguous too;
%     * a sub-structure nested inside a sub-structure, or a struct array.
%
%     * a TOP-LEVEL field whose name belongs to csi_settings or
%       ldpc_settings but which the sim_funs do NOT read from the top level.
%       This is the dangerous case: flattening would give it the same key
%       as the sub-structure field the sim DOES apply, so frames run with
%       the setting ignored would accumulate into a row for the setting
%       applied. Not hypothetical: sim_lookup holds five rows sweeping a
%       top-level max_iters 1..5 that never reached the estimator.
%
%   IDEMPOTENT: an already-flat struct is returned unchanged, so rows from
%   projects that never nested anything keep their existing param_hash.

% Names the sim_funs honour at TOP level as overrides of the sub-structure
% field (checked 2026-10-12: gamma_p in 16 sim_funs, NGS and the PiTau/PiNu/
% Pmax loop in every DD-RELAX/OMP/SAGE sim_fun; ldpc aliases in
% resolve_ldpc_settings). Adding a top-level hook to a sim_fun means adding
% its name here, and nothing else.
top_level_hooks = ["gamma_p","NGS","PiTau","PiNu","Pmax", ...
                   "ldpc_rate","ldpc_N","ldpc_bp_iters","ldpc_colweight"];
% Every csi_settings field any sim_fun reads (grep of Comm Functions,
% 2026-10-12). Only used to recognise a misplaced top-level field, so an
% omission here weakens the guard but never changes a hash.
csi_names = ["method","estimator","max_iters","gamma_p","gamma_p_ref", ...
    "PiTau","PiNu","Pmax","NGS","genie_csi","decoder_aided", ...
    "decoder_aided_mode","decoder_aided_rounds","decoder_aided_sic_rounds", ...
    "decoder_aided_nosub","apg_phase_fix","sigmaI2_perrow","sigmaI2_scalar", ...
    "ridge","ordSel","ordCrit","nAtom","mSmooth","lagWeight","fbAvg", ...
    "fg_trace","dynamic_offset_var","da_probe","sigma_post2_correct", ...
    "row_guard","reliability_sic","pilot_precancel","cancel_mode", ...
    "app_residual_var","adaptive_pmax","op_rounds","full_grid_weight", ...
    "full_grid_refine","decD","PperLag","impl","algo_ver"];

% Already flat: return it untouched. This is what makes the function
% idempotent, and it has to come BEFORE the misplacement guard below: once
% flattened, csi fields such as max_iters legitimately sit at the top level,
% so the guard can only tell a misplaced field from a flattened one while
% the input still has its sub-structures. (Re-hashing a row's stored, flat
% parameters must reproduce its key.)
if ~any(structfun(@isstruct, parameters))
    flat = parameters;
    return;
end

rules = struct( ...
    'csi_settings',  "", ...
    'ldpc_settings', "ldpc_");
renames = struct( ...
    'ldpc_settings', struct('colWeight', "ldpc_colweight"));

flat = struct();
origin = struct();   % flat name -> where it came from, for error messages

top = fieldnames(parameters);
% Plain top-level fields first, so a collision message names both sources.
for i = 1:numel(top)
    v = parameters.(top{i});
    if ~isstruct(v)
        nm = string(top{i});
        if (ismember(nm, csi_names) || startsWith(nm, "ldpc_")) ...
                && ~ismember(nm, top_level_hooks)
            error("flatten_params:unhookedTopLevel", ...
                ['Top-level "%s" is a csi_settings/ldpc_settings field that ' ...
                 'no sim_fun reads from the top level, so the run would ' ...
                 'ignore it while its key claimed it. Move it into the ' ...
                 'sub-structure, or add a top-level hook to the sim_funs ' ...
                 'and list it in top_level_hooks.'], nm);
        end
        flat.(top{i}) = v;
        origin.(top{i}) = string(top{i});
    end
end

for i = 1:numel(top)
    sub = top{i};
    v = parameters.(sub);
    if ~isstruct(v)
        continue;
    end
    if ~isfield(rules, sub)
        error("flatten_params:noRule", ...
            ['parameters.%s is a sub-structure with no flattening rule. ' ...
             'Add one to flatten_params.m that maps its fields onto the ' ...
             'names the sim_funs read, rather than letting it hash nested.'], sub);
    end
    if ~isscalar(v)
        error("flatten_params:structArray", ...
            "parameters.%s is a %s struct array; only scalar structs can be flattened.", ...
            sub, mat2str(size(v)));
    end

    if sub == "ldpc_settings" && ~isfield(v, 'enable')
        v.enable = true;
    end

    fn = fieldnames(v);
    for k = 1:numel(fn)
        if isstruct(v.(fn{k}))
            error("flatten_params:deepNesting", ...
                "parameters.%s.%s is itself a struct; flatten_params handles one level only.", ...
                sub, fn{k});
        end
        name = rules.(sub) + string(fn{k});
        if isfield(renames, sub) && isfield(renames.(sub), fn{k})
            name = renames.(sub).(fn{k});
        end
        if isfield(flat, name)
            error("flatten_params:collision", ...
                ['"%s" is defined by both %s and %s.%s. They are the same ' ...
                 'setting; define it in exactly one place.'], ...
                name, origin.(name), sub, fn{k});
        end
        flat.(name) = v.(fn{k});
        origin.(name) = sub + "." + fn{k};
    end
end
end
