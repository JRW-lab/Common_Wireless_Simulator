function result = cost_head(settings)
%COST_HEAD  Measure what a profile's points COST, as opposed to how well
%they perform. The GUI calls this instead of sim_head whenever a cost
%statistic is selected.
%
%   result = cost_head(settings)
%
%   settings (from WirelessSimulator.buildSettings):
%     profile_sel       index into saved_profiles()
%     figure_statistic  "Runtime: estimator" or "Runtime: detector"
%                       ("FLOPs: ..." is reserved for the analytic-count
%                       pipeline and is refused until that exists)
%     num_frames        timed frames per point; 0 = render the newest
%                       saved measurement instead of measuring
%     use_parallel      must be false when measuring (see below);
%                       rendering a saved result ignores it
%     vis_type          "figure" or "table"
%     legend_loc        as for sim_head (Y limits are ignored; they are
%                       set for error rates)
%     status_fcn        optional @(msg, frac) for GUI progress
%
% WHY THIS EXISTS. The t_ESTcpufull / t_RXcpufull / t_RXfull metrics that
% the sim_funs write to sim_lookup cannot be used to compare cost, and on
% 2026-10-13 the GUI stopped offering them. Each of these alone is enough to
% sink them:
%   1. MACHINES. Desktop and laptop write the same row; sim_save merges by
%      frame-weighted average, and no row records which CPU produced which
%      frames. The stored time is a blend in whatever ratio they raced to.
%   2. CONTENTION. parfor workers compete for cores. Five identical-config
%      rows in the Profile 7 backup spread 1.71-2.52 s on identical work.
%   3. cputime SUMS THREADS. BLAS-heavy stages run at 4-7x their wall time,
%      by a factor that differs between estimators, so the ratio between
%      two estimators is distorted, not just scaled.
%
% WHAT IT DOES INSTEAD (generalising collectors/bench_profile9_timing.m to
% any profile and to the real simulation code):
%   * ONE process, no parallel pool, wall clock (tic/toc) around each stage.
%   * The REAL sim_fun, reached through run_sim_fun.m exactly as collection
%     reaches it, so what is timed is what was collected. Stage times come
%     from the sim_fun's frame_data.t_ESTwall / t_RXwall, which time only
%     the estimator calls and only the detector calls of each frame (all
%     SIC rounds), so frame generation, LDPC decoding and bookkeeping that
%     every arm shares are not in either number.
%   * BACK TO BACK. One warm-up frame per point (JIT, LDPC construction,
%     caches) is discarded. Then each round runs one frame of EVERY point,
%     in a freshly shuffled order, so slow drift in machine state (thermal
%     throttling, a background update) lands on all arms alike instead of
%     on whichever arm happened to run last.
%   * COMMON RANDOM NUMBERS. Every point in round r is run from the same
%     seed, so the arms see the same channel and noise draws; arm-to-arm
%     differences are not channel luck.
%   * Medians with interquartile ranges, because a stage time has a long
%     right tail (an OS interrupt adds, nothing subtracts).
%
% WHAT IT DOES NOT FIX: another program loading the machine. It records the
% conditions with every result (host, thread count, and how many other
% MATLAB processes were using CPU when it started, which is what a
% collector left running looks like) and warns when one is busy. Idle
% MATLAB sessions and idle pool workers are counted but not warned about.
% Measure on an otherwise idle machine.
%
% RESULTS are machine-specific, so they are NOT written to sim_lookup.
% Each run is saved to Logs/runtime/<profile>_<host>_<stamp>.mat with every
% raw sample, and the render path reads the newest one whose points still
% match the profile.

if ~isfield(settings,'status_fcn') || isempty(settings.status_fcn)
    settings.status_fcn = @(msg, frac) [];   % measure/render print to the console themselves
end
stat = string(settings.figure_statistic);
[kind, stage] = parse_statistic(stat);

if kind == "flops"
    error("cost_head:flopsNotYet", ...
        ['"%s" is reserved for the analytic FLOP-count pipeline, which ' ...
         'has not been built yet.'], stat);
end

cws = fileparts(mfilename('fullpath'));
cd(cws);
add_project_paths(cws);

[all_profiles, profile_names] = saved_profiles();
profile = all_profiles{settings.profile_sel};
pname = string(profile_names(settings.profile_sel));
if isfield(profile,'secondary_var') && ~isempty(profile.secondary_var)
    error("cost_head:secondaryVar", ...
        ['Profile "%s" sweeps a secondary variable. The runtime pipeline ' ...
         'handles primary-variable profiles only so far.'], pname);
end

% ---- the profile's points, by the one shared definition
pv = profile.primary_vals;
nP = numel(pv); nC = numel(profile.configs);
pts = struct('ip',{},'ic',{},'par',{},'hash',{});
for ic = 1:nC
    for ip = 1:nP
        [par, h, ok] = build_point_params(profile, pv(ip), ic);
        if ~ok, continue; end
        pts(end+1) = struct('ip',ip,'ic',ic,'par',par,'hash',string(h)); %#ok<AGROW>
    end
end

num_frames = 0;
if isfield(settings,'num_frames'), num_frames = settings.num_frames; end

if num_frames <= 0
    result = load_latest(cws, pname, [pts.hash]);
    msg = sprintf('Rendering runtime measured on %s, %s.', result.host, string(result.when));
    fprintf('%s\n  (saved measurement %s; nothing is being timed)\n', msg, result.file);
    settings.status_fcn(msg, 1);
else
    % ---- refusals: the conditions under which a runtime number is not one.
    % Checked only when MEASURING. Rendering a saved result runs nothing, so
    % an open pool or a ticked Parallelization box cannot affect it.
    if isfield(settings,'use_parallel') && settings.use_parallel
        error("cost_head:parallel", ...
            ['Runtime cannot be measured with Parallelization on: parallel ' ...
             'workers compete for the same cores, which is the contamination ' ...
             'this pipeline exists to avoid. Untick Parallelization and run again.']);
    end
    pool = gcp('nocreate');
    if ~isempty(pool)
        error("cost_head:poolOpen", ...
            ['A parallel pool is still open (%d workers), left over from an ' ...
             'earlier parallel run. Idle workers still hold memory and wake ' ...
             'for housekeeping. Close it with delete(gcp(''nocreate'')) and ' ...
             'run again.'], pool.NumWorkers);
    end

    % Every point is timed over the SAME number of frames, so the arms are
    % compared on equal footing. Frame targets that vary by error rate are
    % a collection setting and do not apply here.
    if isfield(settings,'frame_split') && settings.frame_split
        error("cost_head:frameSplit", ...
            ['Runtime is timed over the same number of frames for every ' ...
             'point. Turn off "Set frames by error rate" and set "All points".']);
    end
    result = measure(pts, num_frames, pname, settings.status_fcn);
    result.primary_var  = profile.primary_var;
    result.primary_vals = pv;
    out = fullfile(cws,'Logs','runtime');
    if ~isfolder(out), mkdir(out); end
    fn = fullfile(out, sprintf('%s_%s_%s.mat', regexprep(char(pname),'[^A-Za-z0-9]+','_'), ...
        regexprep(char(result.host),'[^A-Za-z0-9]+','_'), datestr(result.when,'yyyymmdd_HHMMSS'))); %#ok<DATST>
    save(fn, 'result');
    result.file = string(fn);
    fprintf('\nsaved %s\n', fn);
end

print_summary(result, profile, stage);
if isfield(settings,'vis_type') && string(settings.vis_type) == "table"
    show_table(result, profile, stage);
else
    show_figure(result, profile, stage, settings);
end
end

% =====================================================================
function [kind, stage] = parse_statistic(stat)
switch stat
    case "Runtime: estimator", kind = "runtime"; stage = "est";
    case "Runtime: detector",  kind = "runtime"; stage = "rx";
    case "FLOPs: estimator",   kind = "flops";   stage = "est";
    case "FLOPs: receiver",    kind = "flops";   stage = "rx";
    otherwise
        error("cost_head:unknownStatistic", ...
            '"%s" is not a cost statistic. sim_head handles performance metrics.', stat);
end
end

function add_project_paths(cws)
% Same folders sim_head adds, so run_sim_fun resolves the same functions.
sub = {'Meta Functions', fullfile('Common-Wireless-Infrastructure','Meta Functions'), ...
    'Comm Functions', 'Comm Functions/Custom Functions', 'Comm Functions/Generation Functions', ...
    'Comm Functions/OFDM Functions', 'Comm Functions/OTFS Functions', ...
    'Comm Functions/OTFS-DD Functions', 'Comm Functions/ODDM Functions', ...
    'Comm Functions/ODDM Functions/DD-RELAX-paper', 'Comm Functions/TODDM Functions', ...
    'Comm Functions/TX RX Functions', 'Comm Functions/Coding Functions'};
for k = 1:numel(sub), addpath(fullfile(cws, sub{k})); end
end

function result = measure(pts, nTimed, pname, status)
nPts = numel(pts);
[~, host] = system('hostname'); host = strtrim(string(host));
[nBusy, nOther] = busy_matlab();
if nBusy > 0
    warning("cost_head:otherMatlabBusy", ...
        ['%d other MATLAB process(es) on this machine are using CPU. If one ' ...
         'is collecting, these timings are contaminated. Saved with the result.'], nBusy);
end

% Console banner: a runtime measurement occupies the machine, and anything
% else started on it now lands in the numbers.
fprintf('\n================ RUNTIME MEASUREMENT UNDERWAY ================\n');
fprintf('  Profile        : %s\n', pname);
fprintf('  Points         : %d, each timed over %d frame(s) after 1 warm-up frame\n', nPts, nTimed);
fprintf('  Recording      : estimator and detector wall time per frame\n');
fprintf('  Machine        : %s, %d computation threads\n', host, maxNumCompThreads);
if nBusy > 0
    fprintf('  WARNING        : %d other MATLAB process(es) are using CPU right now\n', nBusy);
end
fprintf('  Leave this machine otherwise idle until it finishes.\n');
fprintf('==============================================================\n');

% Warm-up: one discarded frame per point. Also proves every point's sim_fun
% reports wall-clock stage times before any real time is spent.
fprintf('Warm-up: %d frame(s), discarded...\n', nPts);
t_warm = tic;
for k = 1:nPts
    status(sprintf('Warm-up %d/%d', k, nPts), 0);
    [~, fd] = run_sim_fun(1, pts(k).par);
    if ~isfield(fd,'t_ESTwall') || ~isfield(fd,'t_RXwall')
        error("cost_head:noWallTimes", ...
            ['The simulation function for config %d of "%s" does not report ' ...
             'wall-clock stage times (frame_data.t_ESTwall / t_RXwall). ' ...
             'Only the SP-LDPC sim_fun has them so far; add the two ' ...
             'timers to that sim_fun the same way before benchmarking it.'], ...
            pts(k).ic, pname);
    end
end
fprintf('Warm-up done in %.0f s. Timing %d round(s) of %d points.\n', toc(t_warm), nTimed, nPts);

est = nan(nPts, nTimed); rx = nan(nPts, nTimed); tot = nan(nPts, nTimed);
order = zeros(nTimed, nPts);
shuffle = RandStream('mt19937ar', 'Seed', 20261013);   % order only, not the frames
seed0 = 1e6;
callerRng = rng;
cleanup = onCleanup(@() rng(callerRng));
t_all = tic;
for r = 1:nTimed
    order(r,:) = randperm(shuffle, nPts);
    for j = 1:nPts
        k = order(r,j);
        rng(seed0 + r);                       % same draws for every arm in round r
        t0 = tic;
        [~, fd] = run_sim_fun(1, pts(k).par);
        tot(k,r) = toc(t0);
        est(k,r) = fd.t_ESTwall(1);
        rx(k,r)  = fd.t_RXwall(1);
    end
    el = toc(t_all);
    msg = sprintf('Round %d/%d (%d points each) - about %.0f s left', ...
        r, nTimed, nPts, el/r*(nTimed-r));
    fprintf('  %s\n', msg);
    status(msg, r/nTimed);
end
fprintf('Runtime measurement finished in %.0f s.\n', toc(t_all));

result = struct();
result.profile   = pname;
result.when      = datetime('now');
result.host      = host;
result.threads   = maxNumCompThreads;
result.matlab    = string(version);
result.n_other_matlab = nOther;   % at the start of the run
result.n_busy_matlab  = nBusy;    % of those, using CPU
result.n_timed   = nTimed;
result.n_warmup  = 1;
result.seed0     = seed0;
result.order     = order;
result.ic        = [pts.ic].';
result.ip        = [pts.ip].';
result.hash      = [pts.hash].';
result.est       = est;      % s, per point x round
result.rx        = rx;
result.total     = tot;      % whole sim_fun call, for reference only
end

function [nBusy, nOther] = busy_matlab()
% Other MATLAB processes on this machine, and how many of them are using
% CPU. A bare count is the wrong test: an idle GUI session and its idle
% pool workers cost nothing, while one collecting worker ruins the run.
% Samples two seconds of CPU time per process, excluding this one.
nBusy = NaN; nOther = NaN;
if ~ispc, return; end
try
    me = feature('getpid');
    cmd = ['powershell -NoProfile -Command "$a = Get-Process MATLAB -EA 0 | Select Id, CPU; ' ...
           'Start-Sleep 2; $b = Get-Process MATLAB -EA 0; foreach ($p in $a) { ' ...
           '$q = $b | ? Id -eq $p.Id; if ($q) { ''{0} {1}'' -f $p.Id, ($q.CPU - $p.CPU) } }"'];
    [~, out] = system(cmd);
    v = sscanf(out, '%f %f', [2 Inf]);
    v = v(:, v(1,:) ~= me);
    nOther = size(v, 2);
    nBusy = nnz(v(2,:) > 0.2);    % >10% of one core over the 2 s window
catch
end
end

function result = load_latest(cws, pname, hashes)
d = dir(fullfile(cws,'Logs','runtime', [regexprep(char(pname),'[^A-Za-z0-9]+','_') '_*.mat']));
if isempty(d)
    error("cost_head:noMeasurement", ...
        ['No runtime measurement saved for "%s". Set Number of Trials to the ' ...
         'frames per point to time, untick Parallelization, and press Simulate.'], pname);
end
[~, o] = sort([d.datenum], 'descend');
for k = o
    S = load(fullfile(d(k).folder, d(k).name), 'result');
    if isequal(sort(S.result.hash), sort(hashes(:)))
        result = S.result;
        result.file = string(fullfile(d(k).folder, d(k).name));
        return;
    end
end
error("cost_head:staleMeasurement", ...
    ['Every saved runtime measurement for "%s" was taken on a different set ' ...
     'of points than the profile now defines (it has been edited since). ' ...
     'Run a new measurement.'], pname);
end

function [med, lo, hi] = stats(result, stage)
X = result.(stage) * 1e3;    % ms
med = median(X, 2); lo = prctile_(X, 25); hi = prctile_(X, 75);
end

function q = prctile_(X, p)
% Linear-interpolation percentile along rows, without the Statistics Toolbox.
X = sort(X, 2); n = size(X, 2);
pos = 1 + (n - 1) * p / 100;
a = floor(pos); b = min(a + 1, n); f = pos - a;
q = X(:, a) .* (1 - f) + X(:, b) .* f;
end

function lbl = stage_label(stage)
if stage == "est", lbl = "Channel estimator"; else, lbl = "Detector (PT-MMSE)"; end
end

function print_summary(result, profile, stage)
[med, lo, hi] = stats(result, stage);
fprintf('\n%s wall time per frame, ms (median [IQR]) -- %s, %d threads, %d timed frames/point\n', ...
    stage_label(stage), result.host, result.threads, result.n_timed);
ref = NaN;
for ic = unique(result.ic).'
    lg = string(profile.legend_vec{min(ic, end)});
    k = find(result.ic == ic);
    s = strjoin(arrayfun(@(j) sprintf('%.1f [%.1f-%.1f]', med(j), lo(j), hi(j)), k, 'UniformOutput', false), '   ');
    m = median(med(k));
    if isnan(ref), ref = m; end
    fprintf('  %-40s %s   (%.2fx of first config)\n', lg, s, m / ref);
end
if result.n_busy_matlab > 0
    fprintf('  WARNING: %d other MATLAB process(es) were using CPU when this started.\n', ...
        result.n_busy_matlab);
end
end

function show_figure(result, profile, stage, settings)
[med, lo, hi] = stats(result, stage);
figure('Name', "Runtime - " + result.profile, 'Color', 'w');
hold on; grid on; box on;
x = result.primary_vals;
hs = gobjects(0);
for ic = unique(result.ic).'
    k = find(result.ic == ic);
    [~, s] = sort(result.ip(k)); k = k(s);
    xs = x(result.ip(k)); xs = xs(:);
    sty = char(profile.line_styles{min(ic, end)});
    h = errorbar(xs, med(k), med(k) - lo(k), hi(k) - med(k), sty, ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'CapSize', 6);
    if isfield(profile,'line_colors')
        c = char(profile.line_colors{min(ic, end)});
        set(h, 'Color', c, 'MarkerFaceColor', c);
    end
    hs(end+1) = h; %#ok<AGROW>
end
set(gca, 'YScale', 'log', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
xlabel(string(result.primary_var), 'Interpreter', 'none');
ylabel(stage_label(stage) + " wall time per frame (ms)");
lg = string(profile.legend_vec(1:numel(hs)));
loc = 'best';
if isfield(settings,'legend_loc') && ~isempty(settings.legend_loc), loc = settings.legend_loc; end
legend(hs, lg, 'Interpreter', 'latex', 'Location', loc);
% The GUI's Y limits are left alone on purpose: they are set for error
% rates, and applied to milliseconds they would clip the whole plot.
note = sprintf('%s, %d threads, %d frames/point, median with IQR, %s', ...
    result.host, result.threads, result.n_timed, datestr(result.when, 'yyyy-mm-dd HH:MM')); %#ok<DATST>
if result.n_busy_matlab > 0
    note = [note sprintf(' - %d other MATLAB busy', result.n_busy_matlab)];
end
title({char(result.profile), note}, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 10);
end

function show_table(result, profile, stage)
[med, lo, hi] = stats(result, stage);
x = result.primary_vals;
cfg = strings(numel(result.ic),1);
for j = 1:numel(result.ic), cfg(j) = string(profile.legend_vec{min(result.ic(j), end)}); end
T = table(cfg, x(result.ip).', med, lo, hi, ...
    'VariableNames', {'Config', char(result.primary_var), 'Median_ms', 'P25_ms', 'P75_ms'});
f = uifigure('Name', "Runtime - " + result.profile);
f.Position(3:4) = [760 420];
uitable(f, 'Data', T, 'Position', [10 10 740 400]);
end
