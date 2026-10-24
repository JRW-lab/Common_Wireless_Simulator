function made = local_store_setup()
%LOCAL_STORE_SETUP  Give a fresh clone a results store that needs no setup.
%
%   made = local_store_setup()
%
%   PURPOSE. Results are written to the stores that db_sources.json (the
%   registry, project root) and Results/<COMPUTER>/stores.json name for
%   this computer. A fresh clone has neither, and without them sim_head
%   would try to reach a MySQL server. On such a clone this function
%   writes both files once, registering this computer with one TEXT store
%   at Results/<COMPUTER>/results_text: plain files, no server, no
%   toolbox. launch_gui and sim_head call it before anything is read or
%   written; it does nothing when either file already exists, or when
%   mysql_local.json shows that MySQL was set up the older way.
%
%   OUTPUT
%     made  true when the files were written by this call.
%
%   LIMITS. Both files are gitignored and machine-specific. To use MySQL or
%   SQLite instead, edit them as README.md describes (Results storage).

% Nothing to do when the computer is already configured, either way
root = fileparts(fileparts(mfilename('fullpath')));
made = false;
reg_file = fullfile(root, 'db_sources.json');
if isfile(reg_file) || isfile(fullfile(root, 'mysql_local.json'))
    return
end

% This computer's name and its Results folder
me = char(host_id());
dir_me = fullfile(root, 'Results', me);
if ~isfolder(dir_me), mkdir(dir_me); end

% Its store list: one text primary inside its own Results folder
stores = struct('writer', me, 'stores', {{struct('kind', 'text', 'priority', 1, ...
    'path', ['Results/' me '/results_text'])}});
write_json(fullfile(dir_me, 'stores.json'), stores);

% The registry: this computer as the only source, with a Phase 3 lane block
reg = struct('sources', {{struct('name', 'local', 'role', 'machine', 'writer', me, 'lane_base', 1000)}});
write_json(reg_file, reg);
made = true;
fprintf('First run: results will be stored as text files in %s\n', fullfile(dir_me, 'results_text'));
end

function write_json(f, s)
% Pretty-printed JSON, written in one go.
fid = fopen(f, 'w');
assert(fid > 0, 'local_store_setup:write', 'Cannot write %s.', f);
fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true));
fclose(fid);
end
