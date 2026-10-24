function names = profile_metrics(profile)
%PROFILE_METRICS  The figure statistics at least one config of a profile records.
%
%   names = profile_metrics(profile)
%
%   PURPOSE. The GUI (Common-Wireless-Infrastructure/WirelessSimulator.m)
%   calls this when a profile is selected and lists only these statistics
%   in its Statistic dropdown, so a metric that no config of the profile
%   writes (for example BER_info on an uncoded profile) is never offered.
%   Each config is merged the way sim_head merges it (default_parameters,
%   or flat_defaults for a flat config, then the config's own fields) and
%   passed to config_metrics; the result is the union, in config_metrics'
%   order.
%
%   INPUT
%     profile  one profile struct from saved_profiles().
%
%   OUTPUT
%     names  string row vector of statistic names.
%
%   LIMITS. The primary and secondary sweep variables are not injected:
%   no rule in config_metrics depends on a swept value.

% Union over configs, in first-seen order
names = strings(1, 0);
flat = [];
if isfield(profile, 'flat_configs'), flat = profile.flat_configs; end
for c = 1:numel(profile.configs)

    % The config's base struct: flat_defaults for a flat config that has them, else default_parameters
    base = profile.default_parameters;
    if ismember(c, flat) && isfield(profile, 'flat_defaults') && ~isempty(profile.flat_defaults)
        base = profile.flat_defaults;
    end
    if isfield(profile, 'MUSIC_settings'), base = mergestructs(base, profile.MUSIC_settings); end

    % Config fields override the base, as in build_point_params
    cfg = profile.configs{c};
    f = fieldnames(cfg);
    for k = 1:numel(f), base.(f{k}) = cfg.(f{k}); end
    names = [names config_metrics(base)]; %#ok<AGROW>
end
names = unique(names, 'stable');
end
