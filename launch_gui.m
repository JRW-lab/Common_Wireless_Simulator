% Launches the shared Wireless Simulator GUI (Common-Wireless-Infrastructure/WirelessSimulator.m)
% for this project. Run this instead of opening an .mlapp file.
%
% ---------------------------------------------------------------------
% WHY THIS DOES MORE THAN ONE addpath NOW (2026-09-21).
%
% There are SIX copies of gen_figure.m across the portfolio and they have
% DIVERGED. CWS's own copy lives at
%     Common Wireless Simulator/Common-Wireless-Infrastructure/Meta Functions/
% and is the only one that knows about flat configs (a config measured once
% and drawn as a horizontal line -- see notes/FEATURE_flat_configs.md). The
% top-level C:\MATLAB Projects\Common-Wireless-Infrastructure copy predates
% that, and also predates the timing-metric generalisation and x_log.
%
% This file used to add only the infrastructure FOLDER, not its
% "Meta Functions" subfolder, so which gen_figure won was decided by
% whatever else happened to be on the user's path -- typically the sibling
% project's stale copy, if they had been working there. The symptom is
% nasty and specific: every profile renders fine EXCEPT the one with flat
% configs, because only that one passes empty hash entries and expects the
% renderer to understand them.
%
% It also never added CWS's own "Meta Functions", so build_point_params
% (the single definition of a test point's parameter hash) resolved to
% nothing until sim_head ran its own addpath.
%
% So: add the folders explicitly, put CWS's copies FIRST (addpath prepends),
% and then VERIFY rather than assume. A wrong-copy failure is otherwise
% silent until it produces a confusing error deep inside a render.
% ---------------------------------------------------------------------

cwsRoot = fileparts(mfilename('fullpath'));

addpath(fullfile(cwsRoot, 'Meta Functions'));
addpath(fullfile(cwsRoot, 'Common-Wireless-Infrastructure', 'Meta Functions'));
addpath(fullfile(cwsRoot, 'Common-Wireless-Infrastructure'));

% --- verify the shared functions resolve INSIDE this project ---
mustBeLocal = {'gen_figure','gen_table','mysql_load','jsonencode_sorted', ...
               'build_point_params','saved_profiles','sim_head'};
bad = {};
for k = 1:numel(mustBeLocal)
    p = which(mustBeLocal{k});
    if isempty(p)
        bad{end+1} = sprintf('  %-22s NOT FOUND', mustBeLocal{k}); %#ok<AGROW>
    elseif ~startsWith(lower(p), lower(cwsRoot))
        bad{end+1} = sprintf('  %-22s resolves OUTSIDE this project:\n      %s', ...
            mustBeLocal{k}, p); %#ok<AGROW>
    end
end
if ~isempty(bad)
    error("launch_gui:wrongCopyOnPath", ...
        ['The GUI would run against code from another project, which will\n' ...
         'fail or silently misbehave (the copies have diverged -- only this\n' ...
         'project''s gen_figure understands flat configs):\n\n%s\n\n' ...
         'Fix: restart MATLAB, cd to\n  %s\nand run launch_gui again. If it\n' ...
         'persists, check for a saved pathdef.m pinning another project''s\n' ...
         'folders ahead of this one.'], strjoin(bad, sprintf('\n')), cwsRoot);
end

WirelessSimulator;
