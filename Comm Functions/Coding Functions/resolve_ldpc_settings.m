function ls = resolve_ldpc_settings(parameters, default_N)
% RESOLVE_LDPC_SETTINGS The one definition of a config's LDPC setup.
%
%   ls = resolve_ldpc_settings(parameters)
%   ls = resolve_ldpc_settings(parameters, default_N)
%
%   LDPC parameters live in parameters.ldpc_settings. They are kept out of
%   csi_settings, which describes CHANNEL ESTIMATION only - the two were
%   mixed, and read inconsistently: the superimposed path took them from
%   csi_settings, the guard-isolated paths took most from the top level and
%   ldpc_interleave from csi_settings, and the defaults disagreed between
%   files (turbo_outer_rounds 2 vs 12, turbo_mode hard vs soft).
%
%   Every field is filled, so callers never test isfield again.
%
%   DECODER-AIDED CHANNEL ESTIMATION LIVES IN csi_settings, not here.
%   It was called ldpc_turbo_*, then decoder_aided_* here, and moved to
%   csi_settings on 2026-10-07 once the rename made its nature obvious:
%   the decoder's output never reaches the equaliser as a prior, it only
%   builds a better data-cancellation vector so dd_relax gets a cleaner
%   pilot window. Improving Hhat is its only lever, which is channel
%   estimation and therefore csi_settings' business. It still REQUIRES
%   LDPC, and the sim_fun checks that.
%
%   default_N is the caller's geometric codeword length, used when the
%   config does not set one. It differs per scheme because each frame
%   carries a different number of data bits (512 superimposed, 416 for the
%   1L guard, 328 for the 2L guard), so there is no single right default.
%
%   A legacy ldpc_* field left in csi_settings is an ERROR rather than a
%   silent fallback: reading it would keep the two categories tangled, and
%   ignoring it would run a different code than the config describes while
%   still hashing the field.
%
%   INPUTS:
%     parameters: merged per-point parameter struct
%     default_N:  codeword length to use when ldpc_settings.N is unset

% Function setup
if nargin < 2, default_N = 512; end

% Reject LDPC settings left behind in csi_settings
if isfield(parameters,'csi_settings') && isstruct(parameters.csi_settings)
    stray = fieldnames(parameters.csi_settings);
    stray = stray(startsWith(stray,'ldpc'));
    if ~isempty(stray)
        error("resolve_ldpc_settings:ldpcInCsiSettings", ...
            ['csi_settings carries LDPC field(s): %s. LDPC configuration ' ...
             'moved to parameters.ldpc_settings on 2026-10-07; ' ...
             'csi_settings is for channel estimation only. Move them.'], ...
            strjoin(stray', ', '));
    end
end

% Defaults
ls = struct( ...
    'enable',                   false, ...
    'N',                        default_N, ...
    'rate',                     0.5, ...
    'method',                   'random', ...
    'bp_iters',                 30, ...
    'colWeight',                6, ...
    'interleave',               false);

% Overlay the config
if ~isfield(parameters,'ldpc_settings') || ~isstruct(parameters.ldpc_settings)
    return;
end
cfg = parameters.ldpc_settings;

% Presence of the struct means LDPC is on unless explicitly disabled
ls.enable = true;

known = fieldnames(ls);
given = fieldnames(cfg);
moved = given(startsWith(given,'decoder_aided') | startsWith(given,'turbo'));
if ~isempty(moved)
    error("resolve_ldpc_settings:movedToCsiSettings", ...
        ['ldpc_settings carries decoder-aided field(s): %s. These moved to ' ...
         'csi_settings on 2026-10-07 - the mechanism is decoder-aided ' ...
         'CHANNEL ESTIMATION, not a coding option. Move them.'], ...
        strjoin(moved', ', '));
end

unknown = setdiff(given, known);
if ~isempty(unknown)
    error("resolve_ldpc_settings:unknownField", ...
        ['ldpc_settings has unrecognised field(s): %s. Known fields: %s. ' ...
         'An unrecognised field is hashed but never read, which produces ' ...
         'rows that claim a setting the run never applied.'], ...
        strjoin(unknown', ', '), strjoin(known', ', '));
end
for i = 1:numel(given)
    ls.(given{i}) = cfg.(given{i});
end

% SWEEP HOOKS. build_point_params injects primary_var as a TOP-LEVEL field,
% so a profile cannot sweep anything nested -- the value would be hashed and
% never read, producing a "sweep" of identical configurations. These
% prefixed aliases make the sweepable LDPC parameters reachable:
%     p.primary_var = "ldpc_rate"   sweeps ldpc_settings.rate
%     p.primary_var = "ldpc_N"      sweeps ldpc_settings.N
%     p.primary_var = "ldpc_bp_iters" sweeps ldpc_settings.bp_iters
%     p.primary_var = "ldpc_colweight" sweeps ldpc_settings.colWeight
% Prefixed because a bare "N" is already the Doppler-bin count at top level.
% Defining one in BOTH places is an ERROR, mirroring the NGS guard in
% sim_fun_ODDM_DDRELAXSI_PTMMSE: silently preferring one would produce rows
% claiming a setting the run never applied.
alias = {'ldpc_rate','rate'; 'ldpc_N','N'; 'ldpc_bp_iters','bp_iters'; ...
         'ldpc_colweight','colWeight'};
for i = 1:size(alias,1)
    if isfield(parameters, alias{i,1})
        if any(strcmp(given, alias{i,2}))
            error("resolve_ldpc_settings:sweepCollision", ...
                ['%s is set at top level (for a sweep) AND ldpc_settings.%s ' ...
                 'is set in the config. Define it in exactly one place.'], ...
                alias{i,1}, alias{i,2});
        end
        ls.(alias{i,2}) = parameters.(alias{i,1});
    end
end

% Normalise types
ls.enable     = logical(ls.enable);
ls.interleave = logical(ls.interleave);

end
